// Copyright 2008 Isis Innovation Limited
#include "MapPoint.h"
#include "KeyFrame.h"

void MapPoint::RefreshPixelVectors()
{
    //cout<<"MapPoint::RefreshPixelVectors"<<endl;
    KeyFrame &k = *pPatchSourceKF;

    // Find patch pos in KF camera coords
    // Actually this might not exactly correspond to the patch pos!
    // Treat it as a general point on the plane.
    Vector<3> v3PlanePoint_C = k.se3CfromW * v3WorldPos;
//cout<<k.se3CfromW<<endl;
 //cout<<"v3WorldPos "<<v3WorldPos<<endl;
 //cout<<"v3PlanePoint_C "<<v3PlanePoint_C<<endl;
    // Find the height of this above the plane.
    // Assumes the normal is  pointing toward the camera.
    double dCamHeight = fabs(v3PlanePoint_C * v3Normal_NC);//z值
//cout<<"dCamHeight"<<dCamHeight<<endl;
    double dPixelRate = fabs(v3Center_NC * v3Normal_NC);//
//cout<<"dPixelRate"<<dPixelRate<<endl;
    double dOneRightRate = fabs(v3OneRightFromCenter_NC * v3Normal_NC);
//cout<<"dOneRightRate"<<dOneRightRate<<endl;
    double dOneDownRate = fabs(v3OneDownFromCenter_NC * v3Normal_NC);
//cout<<"dOneDownRate"<<dOneDownRate<<endl;
    // Find projections onto plane
    Vector<3> v3CenterOnPlane_C = v3Center_NC * dCamHeight / dPixelRate;//p(c)    根据比例关系求解当前帧世界点的坐标
//cout<<"v3CenterOnPlane_C"<<v3CenterOnPlane_C<<endl;
    Vector<3> v3OneRightOnPlane_C = v3OneRightFromCenter_NC * dCamHeight / dOneRightRate;//p(r)
//cout<<"v3OneRightOnPlane_C"<<v3OneRightOnPlane_C<<endl;
    Vector<3> v3OneDownOnPlane_C = v3OneDownFromCenter_NC * dCamHeight / dOneDownRate;//p(d)
//cout<<"v3OneDownOnPlane_C"<<v3OneDownOnPlane_C<<endl;

//、、、、、、、、、、、、、、、v3CenterOnPlane_C、v3OneRightOnPlane_C、v3OneDownOnPlane_C的z值是一样的

    // Find differences of these projections in the world frame
v3PixelRight_W = k.se3CfromW.get_rotation().inverse() * (v3OneRightOnPlane_C - v3CenterOnPlane_C);
v3PixelDown_W = k.se3CfromW.get_rotation().inverse() * (v3OneDownOnPlane_C - v3CenterOnPlane_C);

//cout<<"v3PixelRight_W:"<<v3PixelRight_W<<endl;
//cout<<"v3PixelDown_W:"<<v3PixelDown_W<<endl;

//exit(1);
}  
