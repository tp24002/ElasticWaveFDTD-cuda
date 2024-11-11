#pragma once

typedef struct {
  int x;
  int y;
  int z;
} DimI3;

typedef struct {
  double x;
  double y;
  double z;
} DimD3;

typedef struct {
  double rho;
  double K;
  double E;
  double G;
  double nu;
  double ramda;
  double khi;
  double gamma;
  double zeta;
  double eta;   // loss factor
  double omega; // first-order natural angular frequency
} Medium;

typedef struct {
  Medium med;
  DimI3 range;
  DimI3 sp;
} Object;

typedef struct {
  DimI3 Txx;
  DimI3 Tyy;
  DimI3 Tzz;
} SigRan;

typedef struct {
  DimI3 Txy;
  DimI3 Tyz;
  DimI3 Tzx;
} TauRan;

typedef struct {
  DimI3 Vx;
  DimI3 Vy;
  DimI3 Vz;
} VelRan;

typedef struct {
  SigRan sr;
  TauRan tr;
  VelRan vr;
} Range;

typedef struct {
  DimI3 pl1;// min ~ range
  DimI3 pl2;// range ~ max
  double ta;
  double fm;
} Pml;

typedef struct {
  double dx;
  double dy;
  double dz;
  double dt;
} Diff;

typedef struct {
  double ramda;
  double mu;
  double c11;
  double rho;
  double zetaxx;
  double zetaxy;
  double zetaxz;
  double zetayx;
  double zetayy;
  double zetayz;
  double zetazx;
  double zetazy;
  double zetazz;
  double gamma;
  double khi;
  double xi11;
  double zetadx;
  double zetady;
  double zetadz;
} MedArr; // 19

typedef struct {
  double *Txx;
  double *Txxx;
  double *Txxy;
  double *Txxz;
  double *Tyy;
  double *Tyyx;
  double *Tyyy;
  double *Tyyz;
  double *Tzz;
  double *Tzzx;
  double *Tzzy;
  double *Tzzz;
} SigArr; // 12

typedef struct {
  double *Txy;
  double *Txyx;
  double *Txyy;
  double *Tyz;
  double *Tyzy;
  double *Tyzz;
  double *Tzx;
  double *Tzxz;
  double *Tzxx;
} TauArr; // 9

typedef struct {
  double *Vx;
  double *Vxx;
  double *Vxy;
  double *Vxz;
  double *Vy;
  double *Vyx;
  double *Vyy;
  double *Vyz;
  double *Vz;
  double *Vzx;
  double *Vzy;
  double *Vzz;
} VelArr; // 12

typedef struct {
  SigArr sa;
  TauArr ta;
  VelArr va;
} BefAft;

typedef struct {
  double Txx;
  double Tyy;
  double Tzz;
} ImpulseArr; // 3

typedef struct {
  DimI3 in;
  int mode;
  double freq;
} Impulse; 


typedef enum { E_AIR = 0, E_CON, E_M_END } E_Mednum;
typedef enum { E_SINE = 100, E_RCOS } E_KIND_OF_IP;
typedef enum { E_HOLE = 200, E_CLACK } E_KIND_OF_DEFECT;