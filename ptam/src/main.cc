// Copyright 2008 Isis Innovation Limited
// This is the main extry point for PTAM
#include <stdlib.h>
#include <iostream>
#include <gvars3/instances.h>
#include "System.h"
#include <TooN/TooN.h>
using namespace std;
using namespace GVars3;

int main()
{
    cout << "  Welcome to PTAM " << endl;
    cout << "  --------------- " << endl;
    cout << "  Parallel tracking and mapping for Small AR workspaces" << endl;
    cout << "  Copyright (C) Isis Innovation Limited 2008 " << endl;
    cout << endl;
    cout << "  Parsing settings.cfg ...." << endl;
    GUI.LoadFile("../config/settings.cfg");

    GUI.StartParserThread(); // Start parsing of the console input
    atexit(GUI.StopParserThread);

//    TooN::Matrix<3> a;
//    a[0]=TooN::makeVector(1,2,3);
//    a[1]=TooN::makeVector(3,2,1);
//    a[2]=TooN::makeVector(1,0,1);
//    a[2].slice<1,2>() = TooN::makeVector(3,4);
//        cout<<"A="<<a<<endl;
//    TooN::Matrix<2> b = a.slice(0,0,2,2);
//    cout<<"B="<<b<<endl;
//    cout<<"slice::"<<a.slice<0,2,3,1>()<<endl;

    try
    {
        System s;
        s.Run();
    }
    catch(CVD::Exceptions::All e)
    {
        cout << endl;
        cout << "!! Failed to run system; got exception. " << endl;
        cout << "   Exception was: " << endl;
        cout << e.what << endl;
    }
}
