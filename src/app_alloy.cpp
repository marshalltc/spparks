/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator

   
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "app_alloy.h"
#include "solve.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{ZERO,VACANT,OCCUPIED,TOP};
enum{DEPOSITION,NNHOP};
enum{NOSWEEP,RANDOM,RASTER,COLOR,COLOR_STRICT};  // from app_lattice.cpp

#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppAlloy::AppAlloy(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  // these can be changed by model choice, see below

  ninteger = 2;
  ndouble = 0;
  delpropensity = 3;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;
  allow_masking = 0;
  numrandom = 0;
  hopstyle = NNHOP;

  create_arrays();
 
  esites = psites = NULL;
  echeck = pcheck = NULL;
  maxevent = 0;
  events = NULL;
  firstevent = NULL;

 
  hopsite = NULL;
  marklist = NULL;
  mark = NULL;

  allocated = 0;

  // default settings for app-specific commands

  depflag = 0;
  barrierflag = 1;

  // statistics



}

/* ---------------------------------------------------------------------- */

AppAlloy::~AppAlloy()
{

  delete [] esites;
  delete [] psites;
  delete [] echeck;
  delete [] pcheck;
  memory->sfree(events);
  memory->destroy(firstevent);

 
 

  delete [] hopsite;
  delete [] marklist;
  memory->destroy(mark);
}

/* ----------------------------------------------------------------------
   input script commands unique to this app
------------------------------------------------------------------------- */

void AppAlloy::input_app(char *command, int narg, char **arg)
{
  if (sites_exist == 0) {
    char str[128];
    sprintf(str,"Cannot use %s command until sites exist",command);
    error->all(FLERR,str);
  }

  if (!allocated) allocate_data();
  allocated = 1;

   if (strcmp(command,"deposition") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal deposition command");
    if (strcmp(arg[0],"off") == 0) {
      if (narg != 1 ) error->all(FLERR,"Illegal deposition command");
      depflag = 0;
      return;
    }

    if (narg != 8) error->all(FLERR,"Illegal deposition command");
    depflag = 1;
    deprate = atof(arg[0]);
    dir[0] = atof(arg[1]);
    dir[1] = atof(arg[2]);
    dir[2] = atof(arg[3]);
    d0 = atof(arg[4]);
    coordlo = atoi(arg[5]);
    coordhi = atoi(arg[6]);
    particleid = atoi(arg[7]);

    if (deprate < 0.0) error->all(FLERR,"Illegal deposition command");
    if (domain->dimension == 2 && (dir[1] >= 0.0 || dir[2] != 0.0))
      error->all(FLERR,"Illegal deposition command");
    if (domain->dimension == 3 && dir[2] >= 0.0)
      error->all(FLERR,"Illegal deposition command");
    if (d0 < 0.0) error->all(FLERR,"Illegal deposition command");
    if (coordlo < 0 || coordhi > maxneigh || coordlo > coordhi)
      error->all(FLERR,"Illegal deposition command");

    double len = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    dir[0] /= len;
    dir[1] /= len;
    dir[2] /= len;

  } else if (strcmp(command,"bond_energy") == 0) {
    if (narg != 3) error->all(FLERR,"Illegal bond_energy command");
	    
	int i = atoi(arg[0]);
	int j = atoi(arg[1]);

	//if(Q[i][j] > 0) error->all(FLERR, "atomic pair already defined");

	Q[i][j] = Q[j][i] = atof(arg[2]);
      
  }else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppAlloy::grow_app()
{
  lattice = iarray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppAlloy::init_app()
{
  if (depflag && nprocs > 1)
    error->all(FLERR,"Cannot perform deposition in parallel");
  if (depflag && nsector > 1)
    error->all(FLERR,"Cannot perform deposition with multiple sectors");

  if (!allocated) allocate_data();
  allocated = 1;

  dimension = domain->dimension;
  dt_sweep = 1.0/maxneigh;

  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (lattice[i] < VACANT || lattice[i] > TOP) flag = 1;
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ----------------------------------------------------------------------
   setup before each run
------------------------------------------------------------------------- */

void AppAlloy::setup_app()
{
  for (int i = 0; i < nlocal+nghost; i++) echeck[i] = pcheck[i] = 0;

  // clear event list

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppAlloy::site_energy(int i)
{
	int site_species = iarray[1][i];
	double total_energy;

	//sum all the bond energies of coordinated atoms to site i
        for (int j = 0; j < numneigh[i]; j++){
		  int neighbor_species = iarray[1][neighbor[i][j]];
		  if (iarray[0][neighbor[i][j]] == 2) total_energy += Q[site_species][neighbor_species];
	}

	return total_energy;
}


double AppAlloy::site_propensity(int i)
{
  int j,k,m,nsites,ihop,nhop1,eflag;
  double einitial,efinal,edelta,probone,proball;

  clear_events(i);

  if (lattice[i] != OCCUPIED) {
    if (depflag && i == 0) {
      add_event(i,-1,deprate,DEPOSITION);
      return deprate;
    } else return 0.0;
  }

  // nhop1 = 1st neigh hops
  // hopsite = all possible hop sites

  nhop1 = 0;
  for (j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == VACANT) hopsite[nhop1++] = neighbor[i][j];
 
  for (ihop = 0; ihop < nhop1; ihop++) {
    j = hopsite[ihop];

    

  proball = 0.0;
  
  int delta = 1;
   

 //   einitial = site_energy(i);
   
  //  lattice[i] = VACANT;
   // lattice[j] = OCCUPIED;
 
//    efinal = site_energy(j);
    
  //  edelta = efinal - einitial;
 //   
  //  lattice[i] = OCCUPIED;
   // lattice[j] = VACANT;

  if (temperature > 0.0){
       if(edelta < 0) probone = 1.0;
       else    probone = exp(-site_energy(i) * t_inverse);
  }else probone = 0.0;
   
    if (probone > 0.0) {
      eflag = NNHOP;
      add_event(i,j,probone,eflag);
      proball += probone;
    }

}
  // add in single deposition event, stored by site 0

  if (depflag && i == 0) {
    add_event(i,-1,deprate,DEPOSITION);
    proball += deprate;
  }

  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */
void AppAlloy::site_event_rejection(int i, RandomPark *random){}


void AppAlloy::site_event(int i, class RandomPark *random)
{
  int j,m,isite;
  // compare prob to threshhold, break when reach it to select event
  // perform event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
    if (ievent < 0) error->one(FLERR,
			       "Did not reach event propensity threshhold");
  }

  // deposition or hop event
  // for deposition event, find site to deposit on
  // after deposition, reset i and j to that site
  // so propensity around it is updated correctly

  if (events[ievent].style == DEPOSITION) {
    m = find_deposition_site(random);
    if (m < 0) return;
    lattice[m] = OCCUPIED;

  //assign particle id to i2 array to distinguish atomic species
  iarray[1][m] = particleid;



    i = j = m;

  } else {
    

    j = events[ievent].destination;
   
    lattice[i] = VACANT;
    lattice[j] = OCCUPIED;

if(ncoord(j) == 0){
lattice[j] = VACANT;
}else {
  //swap particle id values upon diffusion event

  int buffer1 = iarray[1][i];
  int buffer2 = iarray[1][j];

  iarray[1][i] = buffer2;
  iarray[1][j] = buffer1;
}

  }

  // compute propensity changes for self and swap site and their neighs
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  isite = i2site[i];
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  isite = i2site[j];
  if (isite >= 0) {
    propensity[isite] = site_propensity(j);
    esites[nsites++] = isite;
    echeck[isite] = 1;
  }

  
    nsites += neighbor3(i,&esites[nsites]);
    nsites += neighbor3(j,&esites[nsites]);
  
  solve->update(nsites,esites,propensity);


  for (m = 0; m < nsites; m++) echeck[esites[m]] = 0;
}

/* ----------------------------------------------------------------------
   re-compute propensities out to 2nd neighbors of site I
------------------------------------------------------------------------- */

int AppAlloy::neighbor2(int i, int *sites)
{
  int k,kk,m,mm,isite;

  int nsites = 0;

  for (k = 0; k < numneigh[i]; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      sites[nsites++] = isite;
      echeck[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(mm);
	sites[nsites++] = isite;
	echeck[isite] = 1;
      }
    }
  }

  return nsites;
}

/* ----------------------------------------------------------------------
   re-compute propensities out to 3rd neighbors of site I
------------------------------------------------------------------------- */

int AppAlloy::neighbor3(int i, int *sites)
{
  int k,kk,kkk,m,mm,mmm,isite;

  int nsites = 0;

  for (k = 0; k < numneigh[i]; k++) {
    m = neighbor[i][k];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      sites[nsites++] = isite;
      echeck[isite] = 1;
    }
    for (kk = 0; kk < numneigh[m]; kk++) {
      mm = neighbor[m][kk];
      isite = i2site[mm];
      if (isite >= 0 && echeck[isite] == 0) {
	propensity[isite] = site_propensity(mm);
	sites[nsites++] = isite;
	echeck[isite] = 1;
      }
      for (kkk = 0; kkk < numneigh[mm]; kkk++) {
	mmm = neighbor[mm][kkk];
	isite = i2site[mmm];
	if (isite >= 0 && echeck[isite] == 0) {
	  propensity[isite] = site_propensity(mmm);
	  sites[nsites++] = isite;
	  echeck[isite] = 1;
	}
      }
    }
  }

  return nsites;
}



int AppAlloy::ncoord(int i)
{
  int count = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (lattice[neighbor[i][j]] == OCCUPIED) count++;
  return count;
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppAlloy::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppAlloy::add_event(int i, int destination, 
			      double propensity, int eventflag)
{
  // grow event list and setup free list

  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events = 
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;
  events[freeevent].propensity = propensity;
  events[freeevent].destination = destination;
  events[freeevent].style = eventflag;
  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}


 
/* ----------------------------------------------------------------------
   identify a VACANT site to deposit an atom
   return -1 if could not find a suitable site
------------------------------------------------------------------------- */

int AppAlloy::find_deposition_site(RandomPark *random)
{
  // pick a random position at top of box

  double start[3];
  start[0] = domain->boxxlo + domain->xprd*random->uniform();
  if (dimension == 2) {
    start[1] = domain->boxyhi;
    start[2] = 0.0;
  } else {
    start[1] = domain->boxylo + domain->yprd*random->uniform();
    start[2] = domain->boxzhi;
  }

  // for each vacant site:
  // discard site if neighbor count not between coordlo and coordhi
  // find site whose projected distance is closest to start point

  int i,ncount;
  double dist2start;

  int closesite = -1;
  double closedist = 1.0e20;

  for (i = 0; i < nlocal; i++) {
    if (lattice[i] != VACANT) continue;
    ncount = 0;
    for (int j = 0; j < numneigh[i]; j++)
      if (lattice[neighbor[i][j]] == OCCUPIED) ncount++;
    if (ncount < coordlo || ncount > coordhi) continue;

    if (exceed_limit(i,start,dist2start)) continue;
    if (dist2start < closedist) {
      closedist = dist2start;
      closesite = i;
    }
  }



  return closesite;
}

/* ----------------------------------------------------------------------
   test if site M is within normal distance d0 from incident line
   if so, return 0 and dist2start, else return 1
   site M really becomes a periodic image in XY of M, adjusted via iprd/jprd
   dist2start = dist from site M to starting point of incident line
   dist2start is dist along incident line from start point to
     normal projection point of M
------------------------------------------------------------------------- */

int AppAlloy::exceed_limit(int m, double *start, double &dist2start)
{
  int increment,iprd,jprd;

  iprd = jprd = 0;
  double d0sq = d0*d0;

  double distsq = distsq_to_line(m,start,iprd,jprd,dist2start);
  double newdistsq = distsq_to_line(m,start,iprd-1,jprd,dist2start);
  if (newdistsq < distsq) increment = -1;
  else increment = 1;

  iprd += increment;
  newdistsq = distsq_to_line(m,start,iprd,jprd,dist2start);
  while (newdistsq < distsq) {
    distsq = newdistsq;
    iprd += increment;
    newdistsq = distsq_to_line(m,start,iprd,jprd,dist2start);
  }
  iprd -= increment;

  if (dimension == 3) {
    newdistsq = distsq_to_line(m,start,iprd,jprd-1,dist2start);
    if (newdistsq < distsq) increment = -1;
    else increment = 1;

    jprd += increment;
    newdistsq = distsq_to_line(m,start,iprd,jprd,dist2start);
    while (newdistsq < distsq) {
      distsq = newdistsq;
      jprd += increment;
      newdistsq = distsq_to_line(m,start,iprd,jprd,dist2start);
    }
  }
  jprd -= increment;

  if (distsq > d0sq) return 1;
  distsq = distsq_to_line(m,start,iprd,jprd,dist2start);
  return 0;
}

/* ----------------------------------------------------------------------
   compute normal distsq from site M to incident line of deposition
   site M really becomes a periodic image in XY of M, adjusted via iprd/jprd
   also compute and return dist2start
   dist2start = dist from site M to starting point of incident line
   dist2start is dist along incident line from start point to
     normal projection point of M
------------------------------------------------------------------------- */

double AppAlloy::distsq_to_line(int m, double *start,
				    int iprd, int jprd, double &dist2start)
{
  double delta[3],projection[3],offset[3];

  delta[0] = xyz[m][0] + iprd*domain->xprd - start[0];
  delta[1] = xyz[m][1] + jprd*domain->yprd - start[1];
  delta[2] = xyz[m][2] - start[2];
    
  dist2start = dir[0]*delta[0] + dir[1]*delta[1] + dir[2]*delta[2];
  projection[0] = dist2start*dir[0];
  projection[1] = dist2start*dir[1];
  projection[2] = dist2start*dir[2];
  
  offset[0] = delta[0] - projection[0];
  offset[1] = delta[1] - projection[1];
  offset[2] = delta[2] - projection[2];
  return offset[0]*offset[0] + offset[1]*offset[1] + offset[2]*offset[2];
}

/* ----------------------------------------------------------------------
   allocate data structs that have to wait until sites exist
   so that nlocal,nghost,maxneigh are set
------------------------------------------------------------------------- */

void AppAlloy::allocate_data()
{
  int emax = 1 + maxneigh + maxneigh*maxneigh + maxneigh*maxneigh*maxneigh;
  int pmax = 1 + maxneigh;
  esites = new int[2*emax];
  psites = new int[2*pmax];
  

  echeck = new int[nlocal+nghost];
  pcheck = new int[nlocal+nghost];

  memory->create(firstevent,nlocal,"app:firstevent");


  memory->create(Q, 128, 128, "app:Q");


  hopsite = new int[maxneigh*maxneigh + maxneigh];
  marklist = new int[maxneigh*maxneigh];

  mark = NULL;
  
  if (mark)
    for (int i = 0; i < nlocal+nghost; i++) mark[i] = 0;
}
