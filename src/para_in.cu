#include "../header/struct.h"
#include "../header/init.h"
#include "../header/para_in.h"

#include <stdio.h>

void para_in(Coord *region,Coord *center,Coord *con_st,Coord *con_size,Coord *clack_st,Coord *clack_size,Coord *out,Inpaluse *ip,int *outNum,int *tmax){
    //サイズ
    //可変値
    initCoord(region, 100, 100, 100);
    initCoord(con_st, 3, 3, 3);
    initCoord(con_size, 11, 11, 11);
    //2-14で考える
    initCoord(clack_st, 4, 4, 4);
    initCoord(clack_size, 1, 1, 1);
    ip->freq = 2.0e4;
    ip->mode = E_RCOS;//E_SINE,E_RCOS
    //計算用
    Coord pml;
    initCoord(&pml, 32, 32, 32);
    
    // const Coord halfcon[] = {{(con_size[0]->x - 1) / 2,(con_size[0]->y - 1) / 2,(con_size[0]->z - 1) / 2}};
    Coord halfcon;
    initCoord(&halfcon,(con_size->x - 1) / 2,(con_size->y - 1) / 2,(con_size->z - 1) / 2);
    initCoord(center,(region->x + 2 * pml.x - 1) / 2,(region->y + 2 * pml.y - 1) / 2,(region->z + 2 * pml.z - 1) / 2);
    //座標pml
    ip->in.x = center->x;
    ip->in.y = center->y;
    ip->in.z = center->z + halfcon.z;
    *outNum = 6; 
    // initCoord(&out[0],center->x - halfcon.x + 2, center->y, center->z + halfcon.z);
    // initCoord(&out[1],center->x + halfcon.x - 2, center->y, center->z + halfcon.z);
    // initCoord(&out[2],center->x, center->y - halfcon.y + 2, center->z + halfcon.z);
    // initCoord(&out[3],center->x, center->y + halfcon.y - 2, center->z + halfcon.z);
    initCoord(&out[0],center->x - halfcon.x, center->y, center->z);
    initCoord(&out[1],center->x + halfcon.x, center->y, center->z);
    initCoord(&out[2],center->x, center->y - halfcon.y, center->z);
    initCoord(&out[3],center->x, center->y + halfcon.y, center->z);
    initCoord(&out[4],center->x, center->y, center->z - halfcon.z);
    initCoord(&out[5],center->x, center->y, center->z + halfcon.z);
    // int tmax = 32768;
    *tmax = 16384;
}