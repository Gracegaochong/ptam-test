// Copyright 2008 Isis Innovation Limited
#include "System.h"
#include "OpenGL.h"
#include <gvars3/instances.h>
#include <stdlib.h>
#include "ATANCamera.h"
#include "MapMaker.h"
#include "Tracker.h"
#include "ARDriver.h"
#include "MapViewer.h"
#include <TooN/TooN.h>
#include <vector>
#include <TooN/Cholesky.h>
#include <TooN/wls.h>
#include "HomographyInit.h"
using namespace TooN;
using namespace CVD;
using namespace std;
using namespace GVars3;



System::System()
    : mGLWindow(mVideoSource.Size(), "PTAM")
{
    GUI.RegisterCommand("exit", GUICommandCallBack, this);
    GUI.RegisterCommand("quit", GUICommandCallBack, this);

    mimFrameBW.resize(mVideoSource.Size());
    mimFrameRGB.resize(mVideoSource.Size());
    // First, check if the camera is calibrated.
    // If not, we need to run the calibration widget.
    Vector<NUMTRACKERCAMPARAMETERS> vTest;

    vTest = GV3::get<Vector<NUMTRACKERCAMPARAMETERS> >("Camera.Parameters", ATANCamera::mvDefaultParams, HIDDEN);
    mpCamera = new ATANCamera("Camera");
    Vector<2> v2;
    if(v2==v2) ;
    if(vTest == ATANCamera::mvDefaultParams)
    {
        cout << endl;
        cout << "! Camera.Parameters is not set, need to run the CameraCalibrator tool" << endl;
        cout << "  and/or put the Camera.Parameters= line into the appropriate .cfg file." << endl;
        exit(1);
    }

    mpMap = new Map;
    mpMapMaker = new MapMaker(*mpMap, *mpCamera);
    mpTracker = new Tracker(mVideoSource.Size(), *mpCamera, *mpMap, *mpMapMaker);
    mpARDriver = new ARDriver(*mpCamera, mVideoSource.Size(), mGLWindow);
    mpMapViewer = new MapViewer(*mpMap, mGLWindow);

    GUI.ParseLine("GLWindow.AddMenu Menu Menu");
    GUI.ParseLine("Menu.ShowMenu Root");
    GUI.ParseLine("Menu.AddMenuButton Root Reset Reset Root");
    GUI.ParseLine("Menu.AddMenuButton Root Spacebar PokeTracker Root");
    GUI.ParseLine("DrawAR=0");
    GUI.ParseLine("DrawMap=0");
    GUI.ParseLine("Menu.AddMenuToggle Root \"View Map\" DrawMap Root");
    GUI.ParseLine("Menu.AddMenuToggle Root \"Draw AR\" DrawAR Root");

    mbDone = false;
};

void System::Run()
{

//     Vector<6> trans=makeVector(1,2,3,4,5,6);
//      Vector<3> a=trans.slice<3,3>();
//      cout<<a<<endl;
//      exit(1);
//    Matrix<3> rot(Data(1, 0, 0,
//                      0, 1, 0,
//                      0, 0, 1));
//   Vector<3> trans=makeVector(10,4,-9);
//   SE3<> tform;
//   tform.get_rotation()=rot;
//   tform.get_translation()=trans;
//   Vector<6> mv6SBIRot = tform.ln();//Vector<6> mv6SBIRot 只是得到旋转 还没有得到平移量
//   cout<<mv6SBIRot<<endl;
//    exit(1);
    /*  Vector<3> v3Cam =makeVector(1,2,3);

    Vector<3> v3Motion0 = SO3<>::generator_field(0, v3Cam);
    Vector<3> v3Motion1 = SO3<>::generator_field(1, v3Cam);
    Vector<3> v3Motion2 = SO3<>::generator_field(2, v3Cam);
    cout<<v3Motion0<<endl;
    cout<<v3Motion1<<endl;
    cout<<v3Motion2<<endl;
    Vector<4> v4Cam =makeVector(4,2,3,1);
Vector<4> v4Motion0 = SE3<>::generator_field(0, v4Cam);
Vector<4> v4Motion1 = SE3<>::generator_field(1, v4Cam);
Vector<4> v4Motion2 = SE3<>::generator_field(2, v4Cam);
Vector<4> v4Motion3 = SE3<>::generator_field(3, v4Cam);
Vector<4> v4Motion4 = SE3<>::generator_field(4, v4Cam);
Vector<4> v4Motion5 = SE3<>::generator_field(5, v4Cam);
cout<<endl;
cout<<v4Motion0<<endl;
cout<<v4Motion1<<endl;
cout<<v4Motion2<<endl;
cout<<v4Motion3<<endl;
cout<<v4Motion4<<endl;
cout<<v4Motion5<<endl;
 exit(1);*/
//    SO3<> so3;//3*3
//    SE3<> tform;//3*4
//    exit(1);
//        Matrix<3> rot(Data(0.994197, 0.00196337, -0.107559,
//                           -0.00677516, 0.998991, -0.0443892,
//                           0.107363, 0.0448603, 0.993207));
//        Vector<3> trans=makeVector(10.6,4.10671,-4.25499);
//        SE3<> tform;
//        tform.get_rotation()=rot;
//        tform.get_translation()=trans;
//Vector<3> Wpoint1 = makeVector(-16.5779, -9.41917, 71.9524 );
//Vector<3> v3Cam=tform*Wpoint1;
////Vector<3> v3Cam=makeVector(-59.0747, 82.2991, 17.2184  );
//cout<<"v3Cam"<<v3Cam<<endl;
//Vector<2> v2ImPlane = project(v3Cam);
//Vector<2> v2Image   = mpCamera->Project(v2ImPlane);
//cout << "v2Image: " << v2Image << std::endl;
//Matrix<2,2> m2CamDerivs = mpCamera->GetProjectionDerivs();
//cout << "m2CamDerivs: " << m2CamDerivs << std::endl;
//exit(1);
//     Matrix<3> m3VStar(Data(177, 7.22, 60.6,
//                            7.22, 159, -62.8,
//                            60.6, -62.8, 47.5));
//    Cholesky<3> chol(m3VStar);
//    Matrix<3> m3VStarInv = chol.get_inverse();
//    cout<<"m3VStarInv"<<m3VStarInv<<endl;
//   exit(1);
//    Vector<4> v4Cam = makeVector(1,2,3,1);
//    for(int m=0;m<6;m++)
//    {
//        Vector<4> v4Motion = SE3<>::generator_field(m, v4Cam);
//        cout<<"m:"<<m<<" : "<<v4Motion<<endl;
//    }
//    exit(1);
//    Matrix<3> rot(Data(0.958802, 0.0183186, -0.283486,
//                        -0.0680735, 0.983659, -0.166674 ,
//                        0.2758, 0.179105, 0.944381));
//    Vector<3> trans=makeVector(8.3252, 12.7398,-9.162);
//    SE3<> tform;
//    tform.get_rotation()=rot;
//    tform.get_translation()=trans;
//    for(int m=0;m<3;m++)
//    {
//        Vector<3> v3Motion = tform.get_rotation().get_matrix().T()[m];
//        cout<<"m:"<<m<<" : "<<v3Motion<<endl;
//    }
//     exit(1);
//    vector<pair<ImageRef, ImageRef> > vTrailMatches;
//    pair<ImageRef, ImageRef> p1;
//    p1.first= ImageRef(425,106);
//    p1.second=ImageRef(312.173,105.005);
//    vTrailMatches.push_back(p1);

//    pair<ImageRef, ImageRef> p2;
//    p2.first= ImageRef(215,106);
//    p2.second=ImageRef(58.1032,77.4303);
//    vTrailMatches.push_back(p2);

//    pair<ImageRef, ImageRef> p3;
//    p3.first= ImageRef(215,254);
//    p3.second=ImageRef(36.4955,266.479);
//    vTrailMatches.push_back(p3);

//    pair<ImageRef, ImageRef> p4;
//    p4.first= ImageRef(425,254);
//    p4.second=ImageRef(305.342,268.363);
//    vTrailMatches.push_back(p4);

//    Vector<2> CamPlane01=mpCamera->UnProject(p1.first);
//    Vector<2> CamPlane02=mpCamera->UnProject(p2.first);
//    Vector<2> CamPlane03=mpCamera->UnProject(p3.first);
//    Vector<2> CamPlane04=mpCamera->UnProject(p4.first);

//    Vector<3> Wpoint1 = makeVector(CamPlane01[0]*75, CamPlane01[1]*75,75);
//    Vector<3> Wpoint2 = makeVector(CamPlane02[0]*75, CamPlane02[1]*75,75);
//    Vector<3> Wpoint3 = makeVector(CamPlane03[0]*75, CamPlane03[1]*75,75);
//    Vector<3> Wpoint4 = makeVector(CamPlane04[0]*75, CamPlane04[1]*75,75);
//    cout<<"point:"<<endl<<Wpoint1<<endl<<Wpoint2<<endl<<Wpoint3<<endl<<Wpoint4<<endl;//世界点
//exit(1);
    
//    Vector<2> CamPlane11=mpCamera->UnProject(p1.second);
//    Vector<2> CamPlane12=mpCamera->UnProject(p2.second);
//    Vector<2> CamPlane13=mpCamera->UnProject(p3.second);
//    Vector<2> CamPlane14=mpCamera->UnProject(p4.second);
//140
//    Matrix<3> rot(Data(0.958802, 0.0183186, -0.283486,
//                    -0.0680735, 0.983659, -0.166674 ,
//                    0.2758, 0.179105, 0.944381));
//    Vector<3> trans=makeVector(8.3252, 12.7398,-9.162);
//    SE3<> tform;
//    tform.get_rotation()=rot;
//    tform.get_translation()=trans;
//    cout<<"unproject(Wpoint1):"<<unproject(Wpoint1)<<endl;
//    Vector<4> v=tform * unproject(Wpoint1);
//    cout<<v<<endl;
//    v=v/v[2];
//    Vector<2> vnear=makeVector(v[0],v[1]);
//    cout<<vnear<<endl;
//    Vector<2> vpixel=mpCamera->Project(vnear);
//    cout<<vpixel<<endl;//投影至像素点对比
//    vector<Vector<2> > vv2Verts;
//  exit(1);



//研究主导平面
//    SE3<> aligen=mpMapMaker->CalcPlaneAligner1();//锟斤拷锟斤拷平锟斤拷
//    cout<<"aligen:"<<endl<<aligen<<endl;
//    Vector<3> point= makeVector(0,0,0);
//    cout<<"0 0 0:"<<aligen * point<<endl;
//    cout<<"0 0 01:"<<aligen.inverse() * point<<endl;
//   exit(1);
    while(!mbDone)
    {

        // We use two versions of each video frame:
        // One black and white (for processing by the tracker etc)
        // and one RGB, for drawing.

        // Grab new video frame...
        mVideoSource.GetAndFillFrameBWandRGB(mimFrameBW, mimFrameRGB);//mimFrameBW, mimFrameRGB图片
        static bool bFirstFrame = true;
        if(bFirstFrame)
        {
            mpARDriver->Init();//锟斤拷图
            bFirstFrame = false;
        }

        mGLWindow.SetupViewport();
        mGLWindow.SetupVideoOrtho();
        mGLWindow.SetupVideoRasterPosAndZoom();

        if(!mpMap->IsGood())
            mpARDriver->Reset();

        static gvar3<int> gvnDrawMap("DrawMap", 0, HIDDEN|SILENT);//锟斤拷锟斤拷
        static gvar3<int> gvnDrawAR("DrawAR", 0, HIDDEN|SILENT);//锟斤拷锟斤拷

        bool bDrawMap = mpMap->IsGood() && *gvnDrawMap;
        bool bDrawAR = mpMap->IsGood() && *gvnDrawAR;

        mpTracker->TrackFrame(mimFrameBW, !bDrawAR && !bDrawMap);

        if(bDrawMap)
            mpMapViewer->DrawMap(mpTracker->GetCurrentPose());
        else if(bDrawAR)
            mpARDriver->Render(mimFrameRGB, mpTracker->GetCurrentPose());

        // mGLWindow.GetMousePoseUpdate();
        string sCaption;
        if(bDrawMap)
            sCaption = mpMapViewer->GetMessageForUser();
        else
            sCaption = mpTracker->GetMessageForUser();
        mGLWindow.DrawCaption(sCaption);
        mGLWindow.DrawMenus();
        mGLWindow.swap_buffers();
        mGLWindow.HandlePendingEvents();
    }
}

void System::GUICommandCallBack(void *ptr, string sCommand, string sParams)
{
    if(sCommand=="quit" || sCommand == "exit")
        static_cast<System*>(ptr)->mbDone = true;
}








