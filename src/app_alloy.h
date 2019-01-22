/* ----------------------------------------------------------------------
app_style alloy allows for the simulation of multiple types of atoms
undergoing diffusion and deposition events. It is subject to the same GPL
license as SPPARKS. This application was created by Timothy Marshall 
marshalltc@vcu.edu.
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(alloy,AppAlloy)

#else

#ifndef SPK_APP_ALLOY_H
#define SPK_APP_ALLOY_H

#include "app_lattice.h"

namespace SPPARKS_NS {



class AppAlloy : public AppLattice {
 

 public:
  AppAlloy(class SPPARKS *, int, char **);
  ~AppAlloy();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();
  double site_energy(int);
  double site_propensity(int);
  void site_event(int, class RandomPark *);

//no rKMC but pure virtual function
  void site_event_rejection(int, class RandomPark *);

 private:
  
  int hopstyle;
  int allocated;
  int *esites,*psites;
  int *echeck,*pcheck;


  int dimension;
  int *lattice;

  struct Event {           // one event for an owned site
    double propensity;     // propensity of this event
    int destination;       // local ID of destination site
    int style;             // nearest-neigh hop or deposition
    int next;              // index of next event for this site
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  int particleid;
  int depflag;             // deposition on or off
  double deprate,thetalo,thetahi;
  double d0;
  int coordlo,coordhi;
  double dir[3];

  int barrierflag;          // energy barriers on or off
  double **hbarrier;
  double **Q;


  int *hopsite;             // list of possible hops for one site
  int *mark;                // flagged sites
  int *marklist;            // list of flagged sites


  int neighbor2(int, int *);
  int neighbor3(int, int *);
  int neighbor4(int, int *);

  
  int ncoord(int);
  void clear_events(int);
  void add_event(int, int, double, int);

  int find_deposition_site(class RandomPark *);
  int exceed_limit(int, double *, double &);
  double distsq_to_line(int, double *, int, int, double &);
  void allocate_data();
};

}

#endif
#endif


