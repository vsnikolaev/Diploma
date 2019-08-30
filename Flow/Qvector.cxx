#include "Qvector.h"
#include "particles.h"
#include <math.h>
#include <vector>

Qvector::Qvector() {
	clear();
}

Qvector::~Qvector() {
	clear();
}

void Qvector::clear() {
	Qx = Qy = 0;
    harm = 1;
	for (int i = 0; i < 2; i++) {
		Qxsubev[i] = 0;
		Qysubev[i] = 0;
	}
}

void Qvector::SetHarm(Int_t h){
    harm = h;
}

Int_t Qvector::GetHarm() {
    return harm;
}

//i==1 Qx, i==2, Qy.
Float_t Qvector::GetComponent(int i) {
	if (i == 1) return Qx;
	if (i == 2) return Qy;
	else return -1;
}

//i==1 sub1, i==2, sub2.
Float_t Qvector::GetQxsub(int i) {
	if (i == 1) return Qxsubev[0];
	if (i == 2) return Qxsubev[1];
	else return -1;
}

//i==1 sub1, i==2, sub2.
Float_t Qvector::GetQysub(int i) {
	if (i == 1) return Qysubev[0];
	if (i == 2) return Qysubev[1];
	else return -1;
}

void Qvector::SetQ(Float_t _Qx, Float_t _Qy) {
   Qx=_Qx;
   Qy=_Qy;
}

void Qvector::Recenter(Float_t _corx, Float_t _cory) {
	Qx -= _corx;
	Qy -= _cory;
}

Float_t Qvector::GetEventPlaneAngle() {
	return atan2(Qy, Qx);
}

//true if correct, false if error occurred
void Qvector::SetQSub(Float_t _Qx1, Float_t _Qy1, Float_t _Qx2, Float_t _Qy2) {
   Qxsubev[0]=_Qx1;
   Qysubev[0]=_Qy1;
   Qxsubev[1]=_Qx2;
   Qysubev[1]=_Qy2;
}

void Qvector::RecenterRes(Float_t _corx1, Float_t _cory1, Float_t _corx2, Float_t _cory2) {
	Qxsubev[0] -= _corx1;
	Qxsubev[1] -= _corx2;
	Qysubev[0] -= _cory1;
	Qysubev[1] -= _cory2;
}

Float_t Qvector::GetResolution(Int_t epharm) {
	return cos(harm * (this->GetSubEPAngle(1,epharm) - this->GetSubEPAngle(2,epharm)));
}

//a==1 first sub event, a==2 second...
Float_t Qvector::GetSubEPAngle(Int_t a, Int_t epharm) {
	if (a == 1 || a == 2)
		return 1 / epharm * TMath::ATan2(Qysubev[a - 1], Qxsubev[a - 1]);
	else return 0;
}

//Solo

Float_t Get_res2sub(Float_t a) {
	return sqrt(2*a);
}

Float_t Get_error_res2sub(Float_t a, Float_t a_er) {
	return sqrt(1.0 / (4.0 * a)) * a_er;
}

//input <cos(km(psia-psib)> km=n
Float_t Get_bessel_resolution_ap(Float_t a, Int_t k) {
	Float_t hi = 0;
	Float_t resolution;
	Float_t err = a / 100;
	for (; hi < 3; hi += 0.0001) {
		if (k==1) resolution = 0.626657*hi - 0.09694*pow(hi, 3) + 0.02754*pow(hi, 4) - 0.002283*pow(hi, 5);		//k==1, v1
		if (k==2) resolution = 0.25*pow(hi, 2) - 0.011414*pow(hi, 3) - 0.034726*pow(hi, 4) + 0.006815*pow(hi, 5);	//k==2, v2
		if (k==3) resolution = 1./16. *sqrt(3.14/2.) * pow(hi, 3.)- 3./256. * sqrt(3.14/2.) * pow(hi, 5) + (3 * sqrt(3.14/2.) * pow(hi, 7))/2048.;
        if (abs(resolution - a) < err)
			break;
	}
	if (hi >= 3) return Get_res2sub(a);
	hi = hi*sqrt(2);
	if (k==1) resolution = 0.626657*hi - 0.09694*pow(hi, 3) + 0.02754*pow(hi, 4) - 0.002283*pow(hi, 5);		//k==1, v1
	if (k==2) resolution = 0.25*pow(hi, 2) - 0.011414*pow(hi, 3) - 0.034726*pow(hi, 4) + 0.006815*pow(hi, 5);	//k==2, v2
	if (k==3) resolution = 1./16. *sqrt(3.14/2.) * pow(hi, 3.)- 3./256. * sqrt(3.14/2.) * pow(hi, 5) + (3 * sqrt(3.14/2.) * pow(hi, 7))/2048.;
	resolution = sqrt(resolution);
	return resolution;
}	


Double_t Get_Res(Double_t _chi, Double_t _harm)
{
    Double_t con = TMath::Sqrt(TMath::Pi() / 2) / 2;
    Double_t arg = _chi * _chi / 4.;
    Double_t order1 = (_harm - 1) / 2.;
    Double_t order2 = (_harm + 1) / 2.;
    Double_t res = con * _chi * exp(-arg) * (ROOT::Math::cyl_bessel_i(order1, arg) + ROOT::Math::cyl_bessel_i(order2, arg));
    return res;
}

Double_t Get_Chi(Double_t _res, Double_t _harm, Int_t accuracy) {
    Double_t chi = 2.0;
    Double_t delta = 1.0;
    for (int i = 0; i < accuracy; i++){
        if (Get_Res(chi, _harm) < _res){
            chi = chi + delta;
        }
        else {
            chi = chi - delta;
        }
        delta = delta / 2.;
    }
    return chi;
}