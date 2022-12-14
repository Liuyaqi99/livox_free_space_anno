#ifndef _FREESPACE_HPP
#define _FREESPACE_HPP

#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <Eigen/Core>
#include <pcl/common/transforms.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/common/common.h>
#include <Eigen/Dense>


#include <vector>
#include <math.h>

using namespace std;

#define SELF_CALI_FRAMES 20

//地面栅格图是二维，xy平面；降采样体素图和去噪体素图是三维，xyz
//地面栅格图大小（但用到的其实是下面带了后缀1的），具体N、D、OFF参照后面体素图说明
#define GND_IMG_NX 150
#define GND_IMG_NY 400
#define GND_IMG_DX 0.2
#define GND_IMG_DY 0.2
#define GND_IMG_OFFX 40
#define GND_IMG_OFFY 40

#define GND_IMG_NX1 24
#define GND_IMG_NY1 20
#define GND_IMG_DX1 4
#define GND_IMG_DY1 4
#define GND_IMG_OFFX1 40
#define GND_IMG_OFFY1 40

//N 降采样体素图大小 
#define DN_SAMPLE_IMG_NX 500
#define DN_SAMPLE_IMG_NY 200
#define DN_SAMPLE_IMG_NZ 100
//D 每个降采样体素格代表的实际大小（单位：m）
#define DN_SAMPLE_IMG_DX 0.4
#define DN_SAMPLE_IMG_DY 0.4
#define DN_SAMPLE_IMG_DZ 0.2

// OFF 代表实际原点相对于投影后原点的各轴偏移量（单位：m）；
// 原来的（-50,-40,-10）【最前右下点】变成了新的原点（0,0,0），原来的原点(0,0,0)变成了（50,40,10）
#define DN_SAMPLE_IMG_OFFX 50 // 对应的x轴实际范围是-50m～150m
#define DN_SAMPLE_IMG_OFFY 40 // 对应的y轴实际范围是-40m～40m
#define DN_SAMPLE_IMG_OFFZ 10// 对应的z轴实际范围是-10m～10m

//N 去噪体素图大小 
//（计算方式：实际范围/体素格大小【如z轴 (10-(-10))/DENOISE_IMG_DZ = 100】）
#define DENOISE_IMG_NX 200 
#define DENOISE_IMG_NY 80 
#define DENOISE_IMG_NZ 100
//D 每个去噪体素格代表的实际大小（单位：m）
#define DENOISE_IMG_DX 1 
#define DENOISE_IMG_DY 1 
#define DENOISE_IMG_DZ 0.2

#define FREE_PI 3.14159265



class LivoxFreeSpace
{
    public:
    LivoxFreeSpace();
    // functions
    int GenerateFreeSpace(float* fPoints1, int pointNum, std::vector<float> & free_space);
    void FreeSpace(float* fPoints, int n, float* free_space, int free_space_n);
    void FreeSpaceFilter(float* free_space_small, int n , std::vector<float> & free_space);
    int GroundSegment(int* pLabel,float *fPoints,int pointNum,float fSearchRadius);

    unsigned char *pVImg;

    ~LivoxFreeSpace();

};


#endif
