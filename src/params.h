// Simulation Parameter is stored in a Struct
typedef struct CfdParam{
    double cfl;
    double mach_inf;
    double attack_angle;
    double gamma;
    double p_inf;
} CfdParam;

typedef struct MeshParam{
    int nelem;
    int nnode;
    int niedge;
    int nbedge;
    int nbgroup;
} MeshParam;

typedef struct BoundaryParam{
    char *bound0;
    char *bound1;
    char *bound2;
    char *bound3;
} BoundaryParam;
