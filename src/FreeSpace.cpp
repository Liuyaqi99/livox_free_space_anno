#include "FreeSpace.hpp"

int filter_x[28]={-1,0,1,-3,-2,2,3,-4,4,-4,4,-5,5,-5,5,-5,5,-1,0,1,-3,-2,2,3,-4,4,-4,4};
int filter_y[28]={-5,-5,-5,-4,-4,-4,-4,-3,-3,-2,-2,-1,-1,0,0,1,1,5,5,5,4,4,4,4,3,3,2,2};
int all_x[89]={-1,0,1, 
                -3,-2,-1,0,1,2,3, 
                -4,-3,-2,-1,0,1,2,3,4, 
                -4,-3,-2,-1,0,1,2,3,4, 
                -5,-4,-3,-2,-1,0,1,2,3,4,5,
                -5,-4,-3,-2,-1,0,1,2,3,4,5,
                -5,-4,-3,-2,-1,0,1,2,3,4,5,
                -1,0,1,
                -3,-2-1,0,1,2,3,
                -4,-3,-2,-1,0,1,2,3,4,
                -4,-3,-2,-1,0,1,2,3,4};
int all_y[89]={-5,-5,-5,
                -4,-4,-4,-4,-4,-4,-4,
                -3,-3,-3,-3,-3,-3,-3,-3,-3,
                -2,-2,-2,-2,-2,-2,-2,-2,-2,
                -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                0,0,0,0,0,0,0,0,0,0,0,
                1,1,1,1,1,1,1,1,1,1,1,
                5,5,5,
                4,4,4,4,4,4,4,
                3,3,3,3,3,3,3,3,3,
                2,2,2,2,2,2,2,2,2};

LivoxFreeSpace::LivoxFreeSpace()
{
    this->pVImg=(unsigned char*)calloc(DN_SAMPLE_IMG_NX*DN_SAMPLE_IMG_NY*DN_SAMPLE_IMG_NZ,sizeof(unsigned char));
}
LivoxFreeSpace::~LivoxFreeSpace()
{
    if(this->pVImg!=NULL)
    {
        free(this->pVImg);
    }
}

void LivoxFreeSpace::FreeSpaceFilter(float* free_space_small, int n , std::vector<float> & free_space)
{
    clock_t t0, t1, t2, t3, t4;
    t0 = clock();
    float pixel_size = 0.2, delta_d_in_r = 0.13, delta_r = 0.15; //delta_d_in_r is smaller than pixel_size and delta_r, to make sure all pixels are covered
    Eigen::MatrixXi src = Eigen::MatrixXi::Zero(100/pixel_size, 100/pixel_size), dst = Eigen::MatrixXi::Zero(100/pixel_size, 100/pixel_size);
    
    std::vector<float> delta_t;
    for (float j = 0.0001; j < 50; j += delta_r) // Prepare the delta theta of different radius
    {
        delta_t.push_back(delta_d_in_r/j);
    } 
    for (int i = 0; i < 360; i++)//i为角度
    {
        float r = min(free_space_small[i], free_space_small[(i + 1) % n]);
        r = min(r, free_space_small[(i - 1 + n) % n]);
        r = sqrt(r);//取r为该扇形栅格及左右两个栅格中的距离最小值
        int k = 0;
        for (float j = 0; j < r - 0.5; j += delta_r)//j为范围
        {
            float dt = delta_t[k++];//距离越远，覆盖生成的delta范围越小
            float theta = (i - 180)*FREE_PI/180.0;//栅格角度的弧度制
            for (float t = theta - 0.01; t < theta + 0.01; t+=dt)
            {
                float x = j*cos(t);//极坐标
                float y = j*sin(t);
                int m =int((50.0 - x) / pixel_size);//极坐标投影到图中的坐标
                int nn =int((50.0 - y) / pixel_size);
                src(m, nn) = 1;//表示该点处可行驶
            }
        }
    }

    t1 = clock();
    for (int i = 0; i < 360; i++)
    {
        for (float j = 0; j < 49; j += delta_r)
        {
            float x = j * cos((i - 180)*FREE_PI/180.0);
            float y = j * sin((i - 180)*FREE_PI/180.0);
            int m =int((50.0 - x) / pixel_size);
            int nn =int((50.0 - y) / pixel_size);
            int theta = int(atan2f(y, x) * 180.0 / FREE_PI + 180.0 + 0.5);
            theta = theta % n;
            float r = min(free_space_small[theta], free_space_small[(theta + 1) % n]);
            r = min(r, free_space_small[(theta - 1 + n) % n]);
            if (r > j*j + 1)//所处距离小于该扇形栅格内可达值，也就是统计所有可达范围内的栅格
            {
                int result = 0;
                for (int k = 0; k < 28; k++)  
                {
                    result += src(m + filter_x[k], nn + filter_y[k]);//看它的小邻域（见图）是否都可行驶
                }
                if (result < 28) // check if position (m, nn) is in free space //如果小邻域不是都可以行驶，则跳出
                    break;
                for (int k = 0; k < 89; k++) 
                {
                    dst(m+all_x[k], nn+all_y[k]) = max(1, dst(m+all_x[k], nn+all_y[k]));//该点处大邻域都标记为可以行驶
                }
                dst(m, nn) = 2;//该点可以行驶
            }
        }
    }


    t2 = clock();

    for (int i = 0; i < dst.rows(); i++)
    {
        for (int j = 0; j < dst.cols(); j++)
        {
            if (dst(i, j) > 0)//可行驶
            {
                float x = (100.0 - i*pixel_size) - 50.0;//回归原坐标点
                float y = (100.0 - j*pixel_size) - 50.0;
                free_space.push_back(x);
                free_space.push_back(y);
                free_space.push_back(255);  //此值为强度值，设为255                  
            }
        }
    }
    t3 = clock();
    // printf("filter time: %f, generate map: %f, conv: %f, fs generate: %f\n\n", 1000.0*(t3 - t0) / CLOCKS_PER_SEC,
    //         1000.0*(t1 - t0) / CLOCKS_PER_SEC, 1000.0*(t2 - t1) / CLOCKS_PER_SEC, 1000.0*(t3 - t2) / CLOCKS_PER_SEC);
}


void LivoxFreeSpace::FreeSpace(float* fPoints, int n, float* free_space, int free_space_n)
{
    int thetaId;
    float distance;

    for(int ii=0; ii < free_space_n; ii++)
    {
        free_space[ii] = 2500;//对于每一个角度，设置距离为2500
    }

    for(int pid=0;pid<n;pid++)//对于每一个非地面点
    {
        //如果这个点z值小于3
        if(fPoints[pid*4+2] < 3) // points of high tree, buildings are rejected
        {
            //滤除车体附近的点（后方挂车的点也可以在这一步进行滤除，-2.5 < x < 2.5，-1.2 < y < 1.2）
            if (abs(fPoints[pid*4 + 1]) < 1.2 && abs(fPoints[pid*4]) < 2.5) // reject the near points of robot
                continue;
            distance = fPoints[pid*4]*fPoints[pid*4] + fPoints[pid*4+1]*fPoints[pid*4+1];//在XY平面上的投影距离平方dis
            thetaId = int((atan2f(fPoints[pid*4+1], fPoints[pid*4]) + FREE_PI) * 180.0 / FREE_PI + 0.5);//度数id
            thetaId = thetaId % free_space_n;
            if(free_space[thetaId] > distance && distance > 1)//取在该扇形栅格内，距离车最近的点，为该扇形栅格可达的最远点
            {
                free_space[thetaId] = distance;
            }
        }
    }

}

//生成可行驶区域的主要函数，fPoints1 float数组类型的原始点，pointNum 原始点数量，free_space 可行驶区域点
int LivoxFreeSpace::GenerateFreeSpace(float* fPoints1, int pointNum, std::vector<float> & free_space)
{
    clock_t t0, t1, t2, t3, t4;
    t0 = clock();
    // down sampling 这一部分是进行去噪和降采样
    float *fPoints2=(float*)calloc(pointNum*4,sizeof(float));//存储降采样后的点
    int *idtrans1=(int*)calloc(pointNum,sizeof(int)); //标记每个点投影到体素图里的序号；属于体素图范围的，其值为在体素图里展开后的一维坐标序号，不属于体素图范围的标记为-1
    int *idtransx=(int*)calloc(pointNum,sizeof(int)); //标记每个点在降采样体素图里的x轴坐标
    int *idtransy=(int*)calloc(pointNum,sizeof(int)); //标记每个点在降采样体素图里的y轴坐标
    int *idtrans2=(int*)calloc(pointNum,sizeof(int)); //标记每个降采样后的点投影到降采样体素图里的序号；

    int pntNum = 0;//降采样后的点数

    this->pVImg=(unsigned char*)calloc(DN_SAMPLE_IMG_NX*DN_SAMPLE_IMG_NY*DN_SAMPLE_IMG_NZ,sizeof(unsigned char));
    //将点云投影为NX*NY*NZ的降采样体素图，该体素是否访问过
    std::vector<int> count(DENOISE_IMG_NX*DENOISE_IMG_NY*DENOISE_IMG_NZ, 0);
    //将点云投影为NX*NY*NZ的去噪体素图，计算每个体素内的点云数目

    for(int pid=0;pid<pointNum;pid++)
    {
        //计算点云中每个点投影到去噪体素图中的三维坐标序号
        //（注意：去噪体素图和降采样体素图表示的实际范围一致（不同的只是两者分辨率大小），
        // 所以原点偏移量也一致，因此这里的偏移量OFF采用的是降采样体素图中的，两者偏移量一样）
        int ix=(fPoints1[pid*4]+DN_SAMPLE_IMG_OFFX)/DENOISE_IMG_DX; 
        int iy=(fPoints1[pid*4+1]+DN_SAMPLE_IMG_OFFY)/DENOISE_IMG_DY; 
        int iz=(fPoints1[pid*4+2]+DN_SAMPLE_IMG_OFFZ)/DENOISE_IMG_DZ;
        idtrans1[pid]=-1;
        if((ix>=0)&&(ix<DENOISE_IMG_NX)&&(iy>=0)&&(iy<DENOISE_IMG_NY)&&(iz>=0)&&(iz<DENOISE_IMG_NZ)) 
        {
            //如果该点在去噪体素图范围内
            int idx = iz*DENOISE_IMG_NX*DENOISE_IMG_NY+iy*DENOISE_IMG_NX+ix;//计算点云投影到去噪体素图中展开后的一维坐标序号
            idtrans1[pid]=idx;//idtrans1标记为该点在去噪体素图里的一维序号
            count[idx]++;//该点在去噪体素图中对应体素内的点云数目++
        }
    }

    for(int pid=0;pid<pointNum;pid++)
    {
        if(idtrans1[pid] > -1 && count[idtrans1[pid]] < 3)
        {
            //如果该点在去噪体素图范围内，但是该点所在体素内的点数小于3，则设该点为0，进行去噪滤波
            fPoints1[pid*4] = 0;
            fPoints1[pid*4 + 1] = 0;
            fPoints1[pid*4 + 2] = 0;

        }
    }

    for(int pid=0;pid<pointNum;pid++)
    {
        //计算点云中每个点投影到降采样体素图中的三维坐标序号
        int ix=(fPoints1[pid*4]+DN_SAMPLE_IMG_OFFX)/DN_SAMPLE_IMG_DX;
        int iy=(fPoints1[pid*4+1]+DN_SAMPLE_IMG_OFFY)/DN_SAMPLE_IMG_DY;
        int iz=(fPoints1[pid*4+2]+DN_SAMPLE_IMG_OFFZ)/DN_SAMPLE_IMG_DZ;

        idtrans1[pid]=-1;
        if((ix>=0)&&(ix<DN_SAMPLE_IMG_NX)&&(iy>=0)&&(iy<DN_SAMPLE_IMG_NY)&&(iz>=0)&&(iz<DN_SAMPLE_IMG_NZ))
        {
            //如果该点在降采样体素图范围内
            //idtrans1标记为该点在降采样体素图里的一维序号
            //idtransx标记为该点在降采样体素图里的x轴坐标，idtransy标记为该点在降采样体素图里的y轴坐标
            idtrans1[pid] = iz*DN_SAMPLE_IMG_NX*DN_SAMPLE_IMG_NY+iy*DN_SAMPLE_IMG_NX+ix;
            idtransx[pid] = ix;
            idtransy[pid] = iy;
            if(pVImg[idtrans1[pid]]==0)//后面这句是官方的注释：没有访问过，肯定栅格内会有重复的，所以fPoints2只取第一个
            {
                //这个体素没有访问过的话，将该点作为代表这个体素的点加入到降采样后的点云中去
                fPoints2[pntNum*4] = fPoints1[pid*4];
                fPoints2[pntNum*4 + 1] = fPoints1[pid*4+1];
                fPoints2[pntNum*4 + 2] = fPoints1[pid*4+2];
                fPoints2[pntNum*4 + 3] = fPoints1[pid*4+3];
                idtrans2[pntNum] = pid;//idtrans2表示降采样后的点在降采样体素图里的一维序号
                pntNum++;//降采样后的点数++
            }
                        
            //如果该点在降采样体素图范围内，标记该点所在的体素为访问过
            pVImg[idtrans1[pid]] = 1;
        }
    }

    t1 = clock();

    int *pLabelGnd=(int*)calloc(pntNum,sizeof(int));
    int ground_point_num = GroundSegment(pLabelGnd, fPoints2, pntNum, 1.0);
    //把降采样后的点fPoints2传进来进行地面分割，pntNum为降采样后的点数，pLabelGnd为是否为地面点的标记，1为地面点，0为非地面点；返回值是地面点的数目

    t2 = clock();

    int agnum = pntNum - ground_point_num;//不是地面点的那些点的数目
    float *fPoints3 = (float*)calloc(agnum*4,sizeof(float));//存储非地面点
    int agcnt=0;
    for(int ii=0;ii<pntNum;ii++)
    {
        if(pLabelGnd[ii]==0)
        {
            fPoints3[agcnt*4]=fPoints2[ii*4];
            fPoints3[agcnt*4+1]=fPoints2[ii*4+1];
            fPoints3[agcnt*4+2]=fPoints2[ii*4+2];
            fPoints3[agcnt*4+3]=fPoints2[ii*4+3];
            agcnt++;
        }
        
    }
    float *free_space_small = (float*)calloc(360,sizeof(float));
    this->FreeSpace(fPoints3, agnum, free_space_small, 360);//传入非地面点的点，非地面点的数目和可达距离图角度分辨率，生成可行驶区域
    this->FreeSpaceFilter(free_space_small, 360, free_space);//传入栅格化的最远可达距离，可达距离图角度分辨率，生成以极坐标形式排列的可行驶区域点云

    free(fPoints2);
    free(idtrans1);
    free(idtrans2);
    free(idtransx);
    free(idtransy);
    free(fPoints3);
    free(pLabelGnd);
    free(this->pVImg);
    free(free_space_small);
    std::vector<int>().swap(count);
    t3 = clock();
    // printf("FreeSpace total Time: %f, Downsample: %f, Ground Segment: %f, FreeSpace: %f\n\n", 1000.0*(t3 - t0) / CLOCKS_PER_SEC, 
    //         1000.0*(t1 - t0) / CLOCKS_PER_SEC, 1000.0*(t2 - t1) / CLOCKS_PER_SEC, 1000.0*(t3 - t2) / CLOCKS_PER_SEC);
}

/*
int LivoxFreeSpace::GroundSegment(int* pLabel,float *fPoints,int pointNum,float fSearchRadius)
Fast ground segmentation using rule-based & plane fitting method 
*/
int LivoxFreeSpace::GroundSegment(int* pLabel,float *fPoints,int pointNum,float fSearchRadius)
{
    int gnum=0;

    float *pGndImg1 = (float*)calloc(GND_IMG_NX1*GND_IMG_NY1,sizeof(float));//地面栅格图，每个栅格存储投影到该栅格中的最低点的z值
    int *tmpLabel1 = (int*)calloc(pointNum,sizeof(int));//标记每个点投影到地面栅格图里的序号；属于地面栅格图范围的，其值为在地面栅格图里展开后的一维坐标序号，不属于图范围的标记为-1
    
    for(int ii=0;ii<GND_IMG_NX1*GND_IMG_NY1;ii++)
    {
        pGndImg1[ii]=100;//将地面栅格图每个栅格赋值为100（以方便后面取每个栅格内的最小值）
    }
    for(int pid=0;pid<pointNum;pid++)
    {
        //计算点云中每个点投影到地面栅格图中的二维坐标序号
        int ix= (fPoints[pid*4]+GND_IMG_OFFX1)/(GND_IMG_DX1+0.000001);
        int iy= (fPoints[pid*4+1]+GND_IMG_OFFY1)/(GND_IMG_DY1+0.000001);
        if(ix<0 || ix>=GND_IMG_NX1 || iy<0 || iy>=GND_IMG_NY1)
        {
            //如果这个点超出了地面分割的范围，则设置label为-1，continue
            tmpLabel1[pid]=-1;
            continue;
        }

        int iid=ix+iy*GND_IMG_NX1;//计算该点在地面栅格图中的一维坐标序号
        tmpLabel1[pid]=iid;//标记label为该点在地面栅格图中的一维坐标序号

        if(pGndImg1[iid]>fPoints[pid*4+2])
        {
            //在每个栅格内，取投影到这个栅格的点云z轴最小值来代表该栅格（也就是投影在该栅格里的最低点）
            pGndImg1[iid]=fPoints[pid*4+2];
        }

    }

    int pnum=0;
    for(int pid=0;pid<pointNum;pid++)
    {
        if(tmpLabel1[pid]>=0)//如果这个点在地面栅格图的范围内
        {
            if(pGndImg1[tmpLabel1[pid]]+0.4>fPoints[pid*4+2])//如果这个点和投影在该栅格里的z轴最小值点(最低点)的z轴距离小于0.4m
            {
                pLabel[pid]=1;//记录为可能为地面的点，标记为1
                pnum++;//可能为地面的点数++
            }
        }
    }
    free(pGndImg1);
    free(tmpLabel1);


    for(int pid=0;pid<pointNum;pid++)
    {
        if(pLabel[pid]==1)//如果是可能为地面的点
        {
            if(fPoints[pid*4+2]>1)//如果该点z值超过1m，则不可能为地面点（此举为了排除较高的平台）
            {
                pLabel[pid]=0;//非地面点标为0
            }
            else if(fPoints[pid*4]*fPoints[pid*4]+fPoints[pid*4+1]*fPoints[pid*4+1]<100)
            {
                //如果该点在xy平面上距离原点小于10m（dis^2<100），且z值超过0.5m，则不可能为地面点
                if(fPoints[pid*4+2]>0.5)
                {
                    pLabel[pid]=0;
                }
            }

        }
        else//如果不是可能为地面的点
        {
            if(fPoints[pid*4]*fPoints[pid*4]+fPoints[pid*4+1]*fPoints[pid*4+1]<400)
            {
                //如果该点在xy平面上距离原点小于20m（dis^2<400），且z值小于0.2m，则可能为地面点
                if(fPoints[pid*4+2]<0.2)
                {
                    pLabel[pid]=1;
                }
            }
        }

    }

    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZI>());//pcl库里的点云,存储地面点
    int in_zone[pointNum] = {0};
    for(int pid=0;pid<pointNum;pid++)
    {
        if(fPoints[pid*4]*fPoints[pid*4]+fPoints[pid*4+1]*fPoints[pid*4+1]<400)
        {
            //如果该点在xy平面上距离原点小于20m（dis^2<400），且可能为地面点，则添加到cloud中去
            in_zone[pid] = 1;//表示该点可能为地面点，已被添加
            if (pLabel[pid]==1)
            {
                pcl::PointXYZI p;
                p.x = fPoints[pid*4];
                p.y = fPoints[pid*4 + 1];
                p.z = fPoints[pid*4 + 2];
                p.intensity = fPoints[pid*4 + 3];
                cloud->points.push_back(p);
            }
        }
    }

    Eigen::Matrix3f cov;//协方差矩阵
	Eigen::Vector4f pc_mean;//质心
    Eigen::MatrixXf normal_;//法向量

	pcl::computeMeanAndCovarianceMatrix(*cloud, cov, pc_mean);//计算地面点云的标准化(normalized)协方差矩阵
	Eigen::JacobiSVD<Eigen::MatrixXf> svd(cov,Eigen::DecompositionOptions::ComputeFullU);//奇异值分解，计算协方差矩阵的特征向量
	normal_ = (svd.matrixU().col(2));//取特征值最小的特征向量作为求解的平面法向量
	Eigen::Vector3f seeds_mean = pc_mean.head<3>();
    //默认平面过原点？
	//  normal.T * [x,y,z] = -d
	float d_ = -(normal_.transpose()*seeds_mean)(0,0); //地面点云质心到平面的距离负值
	float th_dist_d_ = 0.3 - d_; //点到平面距离阈值
    Eigen::MatrixXf points(pointNum, 3);
    for(int k = 0; k < pointNum; k++)//将点云转为Eigen模式
    {
        points.row(k) << fPoints[k*4], fPoints[k*4+ 1], fPoints[k*4+ 2];
    }

    // ground plane model
    Eigen::VectorXf result = points * normal_;//计算每个点到平面距离

    for (int k = 0; k < pointNum; k++)
    {
        if (!in_zone[k])//不是可能的地面点，直接跳
            continue;
        if (result[k] < th_dist_d_) //该点到平面距离和质心到平面距离相差小于0.3,为地面点云
        {
            pLabel[k] = 1;
        }
        else
        {
            pLabel[k] = 0;
        }
        
    }


    gnum=0;//计算地面点云数目
    for(int pid=0;pid<pointNum;pid++)
    {
        if(pLabel[pid]==1)
        {
            gnum++;
        }
    }

    return gnum;
}

