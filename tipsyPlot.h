/*
  ## ##      ########  ######## ######## #### ##    ## ########  ######
  ## ##      ##     ## ##       ##        ##  ###   ## ##       ##    ##
#########    ##     ## ##       ##        ##  ####  ## ##       ##
  ## ##      ##     ## ######   ######    ##  ## ## ## ######    ######
#########    ##     ## ##       ##        ##  ##  #### ##             ##
  ## ##      ##     ## ##       ##        ##  ##   ### ##       ##    ##
  ## ##      ########  ######## ##       #### ##    ## ########  ######
*/
#define TYPE_HEADER 1
#define TYPE_GAS 2
#define TYPE_DARK 3
#define TYPE_STAR 4
#define TYPE_BIN 5
#define VAL_NaN 3.402823466e+38

#define AXIS_X 0
#define AXIS_Y 1
#define AXIS_Z 2

#define ERR_MALLOC_FAIL 1
#define ERR_FILE_OPEN 2
#define ERR_NO_PARTICLES 3
#define ERR_MISSING_ARGS 4
#define ERR_UNKNOWN_PARTICLE 5
#define ERR_UNKNOWN_POINTER_SIZE 6
#define ERR_INVALID_ATTRIBUTE 7

#define WARN_REALLOC_SHRINK 1
#define WARN_REALLOC_DATA_LOSS 2

#define ULINE_START \e[4m
#define ULINE_END \e[0m


/*
 ######  ######## ########  ##     ##  ######  ########  ######
##    ##    ##    ##     ## ##     ## ##    ##    ##    ##    ##
##          ##    ##     ## ##     ## ##          ##    ##
 ######     ##    ########  ##     ## ##          ##     ######
      ##    ##    ##   ##   ##     ## ##          ##          ##
##    ##    ##    ##    ##  ##     ## ##    ##    ##    ##    ##
 ######     ##    ##     ##  #######   ######     ##     ######
*/
typedef struct {
    double simtime;
    int nbodies;
    int ndim;
    int nsph;
    int ndark;
    int nstar;
    int pad;
} header;
typedef struct {
    float mass;
    float pos[3];
    float vel[3];
    float rho;
    float temp;
    float eps;
    float metals;
    float phi;
} gas_particle;
typedef struct {
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    float phi;
} dark_particle;
typedef struct {
    float mass;
    float pos[3];
    float vel[3];
    float metals;
    float tform;
    float eps;
    float phi;
} star_particle;
typedef struct {
    float xmin;
    float xmax;
    float ymin;
    float ymax;
    float zmin;
    float zmax;
    int nloadedsph;     // number of particles currently loaded (not VAL_NaN)
    int nloadeddark;
    int nloadedstar;
} attributes;
typedef struct {
    header* head;
    gas_particle* gas;
    dark_particle* dark;
    star_particle* star;
    attributes* attr;
} tipsy;


// function pointers
typedef float (*flop)(float val1, float val2);
typedef float (*calc_bin)(tipsy* tipsyIn, int type, int p);
//typedef float (*calc_var)(bin_particle* binp);
typedef float (*calc_var)(void* particle, int type);
typedef int (*crop)(void* particle, int type);

typedef struct {
    float xmin;
    float xmax;
    float xval;
    int ngas;
    int ndark;
    int nstar;
    gas_particle gas;
    dark_particle dark;
    star_particle star;
} bin_particle;
typedef struct {
    tipsy* sim;         // pointer to the original tipsy the profile was created from
    calc_bin eqbin;
    int nbins;
    float binwidth;
    bin_particle* bin;
} profile;

typedef struct {
    const char* title;
    const char* label;
    const char* shortname;
    calc_var equation;
    int type;
    float ymin;
    float ymax;
    int nbins;
    int npoints;
    double* profile_xs;
    double* profile_ys;
    double* points_xs;
    double* points_ys;
} derivedvar;

/*
######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##  ######
##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ## ##    ##
##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ## ##
######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##  ######
##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####       ##
##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ### ##    ##
##        #######  ##    ##  ######     ##    ####  #######  ##    ##  ######
*/
int cropThinSliceZ(void* particle, int type);

void hackyTipsySplit(tipsy* tipsyIn, tipsy* tipsyLeft, tipsy* tipsyRight);

float calc_rho(void* particle, int type);
float calc_temp(void* particle, int type);
float calc_pressure(void* particle, int type);
float calc_velx(void* particle, int type);
float calc_entropy(void* particle, int type);
float calc_y(void* particle, int type);
float xpos(tipsy* tipsyIn, int type, int p);

// tipsySimEdit.c
void tipsyCenter(tipsy* tipsyIn);
void tipsyCrop(tipsy* tipsyIn, crop op);
void tipsyScaleShrink(tipsy* tipsyIn, const int xShrink, const int yShrink, const int zShrink);
void tipsyScaleExpand(tipsy* tipsyIn, const float xExpand, const float yExpand, const float zExpand);
void tipsyTesselate(tipsy* tipsyIn, const int xTile, const int yTile, const int zTile);
void tipsyTranslate(tipsy* tipsyIn, const float xShift, const float yShift, const float zShift);

// tipsyMemManage.c
tipsy* tipsyCreate(const double simtime, const int nsph, const int ndark, const int nstar);
void tipsyDestroy(tipsy* tipsyIn);
tipsy* tipsyClone(tipsy* tipsyIn);
void tipsyExtend(tipsy* tipsyIn, const int nNewSPH, const int nNewDark, const int nNewStar);
tipsy* tipsyJoin(tipsy* tipsy1, tipsy* tipsy2);

// tipsyProfile.c
profile* profileCreateParticleSpacing(tipsy* tipsyIn, const int nbinssample, const float xmin, const float xmax, calc_bin xs);
profile* profileCreate(tipsy* tipsyIn, const int nbins, const float min, const float max, calc_bin xs);
void profileDestroy(profile* profileIn);
void initializeDerivedVar(derivedvar* variable, const char label[], const char title[], const char shortname[], calc_var equation, int type);
void calculateDerivedVar(derivedvar* variable, profile* profileIn, int type);
void calculateDerivedVarPoints(derivedvar* variable, profile* profileIn, int type);


// tipsyFileIO.c
tipsy* readTipsyStd(const char filename[]);
float* readTipsyDumpScalar(const char filename[]);
float* readTipsyDumpVector(const char filename[], int axis);
float* readTipsyDumpVectorMagnitude(const char filename[]);
int writeTipsyStd(const char filename[], tipsy* tipsyOut);

// gadgetFileIO.c
tipsy* readGadgetToTipsy(const char filename[]);
int writeGadgetFromTipsy(const char filename[], tipsy* tipsyOut);

// tipsyUtils.c
void autoFindBounds(tipsy* tipsyIn);
void tipsySetDefaults(tipsy* tipsyIn);

// particleFlops.c
void pFlop(void* dest, void* src1, void* src2, int type, flop op);
void pFlopGas(gas_particle* dest, gas_particle* src1, gas_particle* src2, flop op);
void pFlopDark(dark_particle* dest, dark_particle* src1, dark_particle* src2, flop op);
void pFlopStar(star_particle* dest, star_particle* src1, star_particle* src2, flop op);

void vFlopGas(gas_particle* dest, gas_particle* src1, float src2, flop op);
void vFlopDark(dark_particle* dest, dark_particle* src1, float src2, flop op);
void vFlopStar(star_particle* dest, star_particle* src1, float src2, flop op);

float flopCopy(float val1, float val2);
float flopSetZero(float val1, float val2);
float flopAdd(float val1, float val2);
float flopDivide(float val1, float val2);

// tipsyMisc.c
void errorCase(const int errorCode);
void warnCase(const int warningCode);

void printGas(gas_particle* p);
void printHeader(header* h);
void printAttr(attributes* a);

int swapEndianInt(const int valIn);
double swapEndianDouble(const double valIn);
float swapEndianFloat(const float valIn);
void swapEndianBatch(const tipsy* tipsyIn, const int type, const int i);


// misc.c
double findMaxVal(double* arrayIn, int len);
double findMinVal(double* arrayIn, int len);
