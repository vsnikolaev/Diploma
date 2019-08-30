#pragma once

#ifndef MYQVEC
#define MYQVEC

#include "particles.h"
#include <fstream>
#include <iostream>
#include <vector>

class Qvector {
private:
    Int_t harm;
	Float_t Qx;
	Float_t Qy;
    Float_t Qxsubev[2];
	Float_t Qysubev[2];
public:
	Qvector();
	virtual ~Qvector();
	virtual void clear();
    virtual void SetHarm(Int_t h);
    virtual Int_t GetHarm();
	virtual Float_t GetComponent(Int_t i);
    virtual Float_t GetQxsub(Int_t i);
    virtual Float_t GetQysub(Int_t i);
	virtual void SetQ(Float_t _Qx, Float_t _Qy);
	virtual void Recenter(Float_t _corx, Float_t _cory);
	virtual Float_t GetEventPlaneAngle();
	//2random sub
	virtual void SetQSub(Float_t _Qx1, Float_t _Qy1, Float_t _Qx2, Float_t _Qy2);
	virtual void RecenterRes(Float_t _corx1, Float_t _cory1, Float_t _corx2, Float_t _cory2);
	virtual Float_t GetResolution(Int_t epharm);
	virtual Float_t GetSubEPAngle(Int_t a, Int_t epharm);
};

Float_t Get_res2sub(Float_t a);
Float_t Get_error_res2sub(Float_t a, Float_t a_er);
Float_t Get_bessel_resolution_ap(Float_t a);
Double_t Get_Res(const Double_t _chi, const Double_t _harm = 1);
Double_t Get_Chi(Double_t _res, Double_t _harm = 1, Int_t accuracy = 100);

#endif // MYQVEC
