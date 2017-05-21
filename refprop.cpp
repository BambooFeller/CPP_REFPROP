// Example C program written by:
// Chris Muzny
// N.I.S.T.
// Chemical Science and Technology Laboratory
// Physical and Chemical Properties of Fluids Division
// (303) 497-5549 
// chris.muzny@nist.gov

//This program demonstrates explicitly linking the subroutines available in
// refprop.dll.  In order to link this code refprop1.h 
// must be available in the current directory.  When executing refprop.dll must be in the dll
// search path (current directory and $PATH).

#include"windows.h"
#include<stdio.h>
#include"refprop1.h"


// Some constants...

const long refpropcharlength = 255;
const long filepathlength = 255;
const long lengthofreference = 3;
const long errormessagelength = 255;
const long ncmax = 20;		// Note: ncmax is the max number of components
const long numparams = 72;
const long maxcoefs = 50;


class Substance
{
public:
	Substance();
	~Substance();



	// Now use the functions.

	// Refprop variables that need to be defined
	//
	// nc = Number of components in the mixture
	// x[NumberOfComponentsInMixtures] = Mole fraction of each component
	// ierr =  An integer flag defining an error
	// hf[] = a character array defining the fluids in a mixture
	// hrf[] = a character array denoting the reference state
	// herr[] = a character array for storing a string - Error message
	// hfmix[] a character array defining the path to the mixture file

	double x[ncmax], xliq[ncmax], xvap[ncmax], f[ncmax];

	long i, ierr;
	char hf[refpropcharlength*ncmax], hrf[lengthofreference + 1],
		herr[errormessagelength + 1], hfmix[refpropcharlength + 1];



	void showSubstance()
	{
		double wm, ttp, tnbp, tc, pc, dc, zc, acf, dip, rgas;
		long info_index = 1;
		INFOdll(info_index, wm, ttp, tnbp, tc, pc, dc, zc, acf, dip, rgas);
		printf("WM,ACF,DIP,TTP,TNBP   %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", wm, acf, dip, ttp, tnbp);
		printf("TC,PC,DC,RGAS         %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", tc, pc, dc, rgas);
	}


	void initialize()
	{

		printf("initialize the program.\n");


		// First create a pointer to an instance of the library
		// Then have windows load the library.
		HINSTANCE RefpropdllInstance;
		//This looks only in the current directory for refprop.dll
		RefpropdllInstance = LoadLibrary(L"refprop.dll");

		// Then get pointers into the dll to the actual functions.
		ABFL1dll = (fp_ABFL1dllTYPE)GetProcAddress(RefpropdllInstance, "ABFL1dll");
		ABFL2dll = (fp_ABFL2dllTYPE)GetProcAddress(RefpropdllInstance, "ABFL2dll");
		ACTVYdll = (fp_ACTVYdllTYPE)GetProcAddress(RefpropdllInstance, "ACTVYdll");
		AGdll = (fp_AGdllTYPE)GetProcAddress(RefpropdllInstance, "AGdll");
		CCRITdll = (fp_CCRITdllTYPE)GetProcAddress(RefpropdllInstance, "CCRITdll");
		CP0dll = (fp_CP0dllTYPE)GetProcAddress(RefpropdllInstance, "CP0dll");
		CRITPdll = (fp_CRITPdllTYPE)GetProcAddress(RefpropdllInstance, "CRITPdll");
		CSATKdll = (fp_CSATKdllTYPE)GetProcAddress(RefpropdllInstance, "CSATKdll");
		CV2PKdll = (fp_CV2PKdllTYPE)GetProcAddress(RefpropdllInstance, "CV2PKdll");
		CVCPKdll = (fp_CVCPKdllTYPE)GetProcAddress(RefpropdllInstance, "CVCPKdll");
		CVCPdll = (fp_CVCPdllTYPE)GetProcAddress(RefpropdllInstance, "CVCPdll");
		DBDTdll = (fp_DBDTdllTYPE)GetProcAddress(RefpropdllInstance, "DBDTdll");
		DBFL1dll = (fp_DBFL1dllTYPE)GetProcAddress(RefpropdllInstance, "DBFL1dll");
		DBFL2dll = (fp_DBFL2dllTYPE)GetProcAddress(RefpropdllInstance, "DBFL2dll");
		DDDPdll = (fp_DDDPdllTYPE)GetProcAddress(RefpropdllInstance, "DDDPdll");
		DDDTdll = (fp_DDDTdllTYPE)GetProcAddress(RefpropdllInstance, "DDDTdll");
		DEFLSHdll = (fp_DEFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "DEFLSHdll");
		DHD1dll = (fp_DHD1dllTYPE)GetProcAddress(RefpropdllInstance, "DHD1dll");
		DHFLSHdll = (fp_DHFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "DHFLSHdll");
		DIELECdll = (fp_DIELECdllTYPE)GetProcAddress(RefpropdllInstance, "DIELECdll");
		DOTFILLdll = (fp_DOTFILLdllTYPE)GetProcAddress(RefpropdllInstance, "DOTFILLdll");
		DPDD2dll = (fp_DPDD2dllTYPE)GetProcAddress(RefpropdllInstance, "DPDD2dll");
		DPDDKdll = (fp_DPDDKdllTYPE)GetProcAddress(RefpropdllInstance, "DPDDKdll");
		DPDDdll = (fp_DPDDdllTYPE)GetProcAddress(RefpropdllInstance, "DPDDdll");
		DPDTKdll = (fp_DPDTKdllTYPE)GetProcAddress(RefpropdllInstance, "DPDTKdll");
		DPDTdll = (fp_DPDTdllTYPE)GetProcAddress(RefpropdllInstance, "DPDTdll");
		DPTSATKdll = (fp_DPTSATKdllTYPE)GetProcAddress(RefpropdllInstance, "DPTSATKdll");
		DSFLSHdll = (fp_DSFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "DSFLSHdll");
		ENTHALdll = (fp_ENTHALdllTYPE)GetProcAddress(RefpropdllInstance, "ENTHALdll");
		ENTROdll = (fp_ENTROdllTYPE)GetProcAddress(RefpropdllInstance, "ENTROdll");
		ESFLSHdll = (fp_ESFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "ESFLSHdll");
		FGCTYdll = (fp_FGCTYdllTYPE)GetProcAddress(RefpropdllInstance, "FGCTYdll");
		FPVdll = (fp_FPVdllTYPE)GetProcAddress(RefpropdllInstance, "FPVdll");
		GERG04dll = (fp_GERG04dllTYPE)GetProcAddress(RefpropdllInstance, "GERG04dll");
		GETFIJdll = (fp_GETFIJdllTYPE)GetProcAddress(RefpropdllInstance, "GETFIJdll");
		GETKTVdll = (fp_GETKTVdllTYPE)GetProcAddress(RefpropdllInstance, "GETKTVdll");
		GIBBSdll = (fp_GIBBSdllTYPE)GetProcAddress(RefpropdllInstance, "GIBBSdll");
		HSFLSHdll = (fp_HSFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "HSFLSHdll");
		INFOdll = (fp_INFOdllTYPE)GetProcAddress(RefpropdllInstance, "INFOdll");
		LIMITKdll = (fp_LIMITKdllTYPE)GetProcAddress(RefpropdllInstance, "LIMITKdll");
		LIMITSdll = (fp_LIMITSdllTYPE)GetProcAddress(RefpropdllInstance, "LIMITSdll");
		LIMITXdll = (fp_LIMITXdllTYPE)GetProcAddress(RefpropdllInstance, "LIMITXdll");
		MELTPdll = (fp_MELTPdllTYPE)GetProcAddress(RefpropdllInstance, "MELTPdll");
		MELTTdll = (fp_MELTTdllTYPE)GetProcAddress(RefpropdllInstance, "MELTTdll");
		MLTH2Odll = (fp_MLTH2OdllTYPE)GetProcAddress(RefpropdllInstance, "MLTH2Odll");
		NAMEdll = (fp_NAMEdllTYPE)GetProcAddress(RefpropdllInstance, "NAMEdll");
		PDFL1dll = (fp_PDFL1dllTYPE)GetProcAddress(RefpropdllInstance, "PDFL1dll");
		PDFLSHdll = (fp_PDFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "PDFLSHdll");
		PEFLSHdll = (fp_PEFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "PEFLSHdll");
		PHFL1dll = (fp_PHFL1dllTYPE)GetProcAddress(RefpropdllInstance, "PHFL1dll");
		PHFLSHdll = (fp_PHFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "PHFLSHdll");
		PQFLSHdll = (fp_PQFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "PQFLSHdll");
		PREOSdll = (fp_PREOSdllTYPE)GetProcAddress(RefpropdllInstance, "PREOSdll");
		PRESSdll = (fp_PRESSdllTYPE)GetProcAddress(RefpropdllInstance, "PRESSdll");
		PSFL1dll = (fp_PSFL1dllTYPE)GetProcAddress(RefpropdllInstance, "PSFL1dll");
		PSFLSHdll = (fp_PSFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "PSFLSHdll");
		PUREFLDdll = (fp_PUREFLDdllTYPE)GetProcAddress(RefpropdllInstance, "PUREFLDdll");
		QMASSdll = (fp_QMASSdllTYPE)GetProcAddress(RefpropdllInstance, "QMASSdll");
		QMOLEdll = (fp_QMOLEdllTYPE)GetProcAddress(RefpropdllInstance, "QMOLEdll");
		SATDdll = (fp_SATDdllTYPE)GetProcAddress(RefpropdllInstance, "SATDdll");
		SATEdll = (fp_SATEdllTYPE)GetProcAddress(RefpropdllInstance, "SATEdll");
		SATHdll = (fp_SATHdllTYPE)GetProcAddress(RefpropdllInstance, "SATHdll");
		SATPdll = (fp_SATPdllTYPE)GetProcAddress(RefpropdllInstance, "SATPdll");
		SATSdll = (fp_SATSdllTYPE)GetProcAddress(RefpropdllInstance, "SATSdll");
		SATTdll = (fp_SATTdllTYPE)GetProcAddress(RefpropdllInstance, "SATTdll");
		SETAGAdll = (fp_SETAGAdllTYPE)GetProcAddress(RefpropdllInstance, "SETAGAdll");
		SETKTVdll = (fp_SETKTVdllTYPE)GetProcAddress(RefpropdllInstance, "SETKTVdll");
		SETMIXdll = (fp_SETMIXdllTYPE)GetProcAddress(RefpropdllInstance, "SETMIXdll");
		SETMODdll = (fp_SETMODdllTYPE)GetProcAddress(RefpropdllInstance, "SETMODdll");
		SETREFdll = (fp_SETREFdllTYPE)GetProcAddress(RefpropdllInstance, "SETREFdll");
		SETUPdll = (fp_SETUPdllTYPE)GetProcAddress(RefpropdllInstance, "SETUPdll");
		SPECGRdll = (fp_SPECGRdllTYPE)GetProcAddress(RefpropdllInstance, "SPECGRdll");
		SUBLPdll = (fp_SUBLPdllTYPE)GetProcAddress(RefpropdllInstance, "SUBLPdll");
		SUBLTdll = (fp_SUBLTdllTYPE)GetProcAddress(RefpropdllInstance, "SUBLTdll");
		SURFTdll = (fp_SURFTdllTYPE)GetProcAddress(RefpropdllInstance, "SURFTdll");
		SURTENdll = (fp_SURTENdllTYPE)GetProcAddress(RefpropdllInstance, "SURTENdll");
		TDFLSHdll = (fp_TDFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "TDFLSHdll");
		TEFLSHdll = (fp_TEFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "TEFLSHdll");
		THERM0dll = (fp_THERM0dllTYPE)GetProcAddress(RefpropdllInstance, "THERM0dll");
		THERM2dll = (fp_THERM2dllTYPE)GetProcAddress(RefpropdllInstance, "THERM2dll");
		THERM3dll = (fp_THERM3dllTYPE)GetProcAddress(RefpropdllInstance, "THERM3dll");
		THERMdll = (fp_THERMdllTYPE)GetProcAddress(RefpropdllInstance, "THERMdll");
		THFLSHdll = (fp_THFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "THFLSHdll");
		TPFLSHdll = (fp_TPFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "TPFLSHdll");
		TPRHOdll = (fp_TPRHOdllTYPE)GetProcAddress(RefpropdllInstance, "TPRHOdll");
		TQFLSHdll = (fp_TQFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "TQFLSHdll");
		TRNPRPdll = (fp_TRNPRPdllTYPE)GetProcAddress(RefpropdllInstance, "TRNPRPdll");
		TSFLSHdll = (fp_TSFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "TSFLSHdll");
		VIRBdll = (fp_VIRBdllTYPE)GetProcAddress(RefpropdllInstance, "VIRBdll");
		VIRCdll = (fp_VIRCdllTYPE)GetProcAddress(RefpropdllInstance, "VIRCdll");
		WMOLdll = (fp_WMOLdllTYPE)GetProcAddress(RefpropdllInstance, "WMOLdll");
		XMASSdll = (fp_XMASSdllTYPE)GetProcAddress(RefpropdllInstance, "XMASSdll");
		XMOLEdll = (fp_XMOLEdllTYPE)GetProcAddress(RefpropdllInstance, "XMOLEdll");


		//注释
		//混合物调用使用以下方法，而不是以上方法
		//...For a mixture, use the following setup instead of the lines above.
		// Use "|" as the file name delimiter for mixtures
		//   i=3;
		//   strcpy(hf,"nitrogen.fld");
		//strcat(hf,"|argon.fld");
		//strcat(hf,"|oxygen.fld");
		//   strcpy(hfmix,"hmx.bnc");
		//   strcpy(hrf,"DEF");
		//strcpy(herr,"Ok");
		//x[0]=.7812;     //Air composition
		//   x[1]=.0092;
		//   x[2]=.2096;
		/* ********************************************/



		////...Call SETUP to initialize the program
		//SETUPdll(i, hf, hfmix, hrf, ierr, herr,
		//	refpropcharlength*ncmax, refpropcharlength,
		//	lengthofreference, errormessagelength);
		//if (ierr != 0) printf("%s\n", herr);

		//double wm, ttp, tnbp, tc, pc, dc, zc, acf, dip, rgas;
		//long info_index = 1;
		//INFOdll(info_index, wm, ttp, tnbp, tc, pc, dc, zc, acf, dip, rgas);
		//printf("WM,ACF,DIP,TTP,TNBP   %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", wm, acf, dip, ttp, tnbp);
		//printf("TC,PC,DC,RGAS         %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", tc, pc, dc, rgas);
		////...Calculate molecular weight of a mixture
		////     wm=WMOLdll(x)

		////...Get saturation properties given t,x; for i=1: x is liquid phase
		////.....                                   for i=2: x is vapor phase

		//double t = 100.0;
		//double p, dl, dv;

		//SATTdll(t, x, i, p, dl, dv, xliq, xvap, ierr, herr, errormessagelength);
		//printf("P,Dl,Dv,xl[0],xv[0]   %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", p, dl, dv, xliq[0], xvap[0]);
		//i = 2;
		//SATTdll(t, x, i, p, dl, dv, xliq, xvap, ierr, herr, errormessagelength);
		//printf("P,Dl,Dv,xl[0],xv[0]   %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", p, dl, dv, xliq[0], xvap[0]);

		////...Calculate saturation properties at a given p. i is same as SATT

		//i = 2;
		//SATPdll(p, x, i, t, dl, dv, xliq, xvap, ierr, herr, errormessagelength);
		//printf("T,Dl,Dv,xl(1),xv(1)   %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", t, dl, dv, xliq[0], xvap[0]);

		////...Other saturation routines are given in SAT_SUB.FOR

		//t = 300.0;
		//p = 2000.0;

		//t = 100;
		//p = 100;
		//x[0] = 1;
		//i = 1; //i = 1 ,液体；i = 2,气体
		//printf("给定温度t = %f，压力 p = %f，x = %f \n", t, p, x);
		////...Calculate d from t,p,x
		////...If phase is known: (j=1: Liquid, j=2: Vapor)
		//long j = 2;
		//double d, q, e, h, s, cv, cp, w, b, c,
		//	dpdrho, d2pdd2, dpdt, dhdt_d, dhdt_p, dhdp_t, dhdp_d,
		//	sigma, dhdd_t, dhdd_p, eta, tcx, pp, tt, hjt, h1, dd;
		//long tmp_int = 0;
		//TPRHOdll(t, p, x, j, tmp_int, d, ierr, herr, errormessagelength);
		//printf("T,P,D                 %10.4f,%10.4f,%10.4f\n", t, p, d);

		////...If phase is not known, call TPFLSH
		////...Calls to TPFLSH are much slower than TPRHO since SATT must be called first.
		////.....(If two phase, quality is returned as q)

		////
		//printf("不清楚气态还是液态\n");
		//t = 76;
		//p = 1000;
		//i = 1;
		////x = 1.0;
		//printf("第二种t =%f ,p = %f,x[0] = %f \n", t, p, x[0]);
		//TPFLSHdll(t, p, x, d, dl, dv, xliq, xvap, q, e, h, s, cv, cp, w, ierr, herr, errormessagelength);
		//printf("T,P,D,H,CP            %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", t, p, d, h, cp);
		//printf("密度d= %f (kg/m3)", d*wm);

		////...Calculate pressure (p), internal energy (e), enthalpy (h), entropy (s),
		////.....isochoric (cv) and isobaric (cp) heat capacities, speed of sound (w),
		////.....and Joule-Thomson coefficient (hjt) from t,d,x
		////.....(subroutines THERM2 and THERM3 contain more properties, see PROP_SUB.FOR)
		//THERMdll(t, d, x, p, e, h, s, cv, cp, w, hjt);

		////...Calculate pressure
		//PRESSdll(t, d, x, p);

		////...Calculate fugacity
		//FGCTYdll(t, d, x, f);

		////...Calculate second and third virial coefficients
		//VIRBdll(t, x, b);
		//VIRCdll(t, x, c);
		//printf("F,B,C                 %10.4f,%10.4f,%10.4f\n", f[0], b, c);

		////...Calculate the derivatives: dP/dD, d^2P/dD^2, dP/dT  (D indicates density)
		////...(dD/dP, dD/dT, and dB/dT are also available, see PROP_SUB.FOR)
		//DPDDdll(t, d, x, dpdrho);
		//DPDD2dll(t, d, x, d2pdd2);
		//DPDTdll(t, d, x, dpdt);
		//printf("dP/dD,d2P/dD2,dP/dT   %10.4f,%10.4f,%10.4f\n", dpdrho, d2pdd2, dpdt);


		////...Calculate derivatives of enthalpy with respect to T, P, and D
		//DHD1dll(t, d, x, dhdt_d, dhdt_p, dhdd_t, dhdd_p, dhdp_t, dhdp_d);
		//printf("Enthalpy derivatives  %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n",
		//	dhdt_d, dhdt_p, dhdd_t, dhdd_p / 1000.0, dhdp_t);
		////...Calculate surface tension
		//SURFTdll(t, dl, x, sigma, ierr, herr, errormessagelength);
		//printf("T,SURF. TN.           %10.4f,%10.4f\n", t, sigma);

		////...Calculate viscosity (eta) and thermal conductivity (tcx)
		//printf("以下计算粘度（(10^-6 Pa.s)和导热系数（mW/(m.K)）\n");
		//TRNPRPdll(t, d, x, eta, tcx, ierr, herr, errormessagelength);
		//printf("VIS.,TH.CND.          %10.4f,%10.4f\n", eta, tcx*1000.0);

		//printf("给定 t,d,x 计算物性\n");
		//t = 77;
		//printf("t = %f, d =  %f x = %f\n", t, d, x);
		////...General property calculation with inputs of t,d,x
		//TDFLSHdll(t, d, x, pp, dl, dv, xliq, xvap, q, e, h1, s, cv, cp, w, ierr, herr, errormessagelength);
		//printf("T, D, P from TDFLSH   %10.4f,%10.4f,%10.4f\n", t, d, pp / 1000.0);

		////...General property calculation with inputs of p,d,x
		//PDFLSHdll(p, d, x, tt, dl, dv, xliq, xvap, q, e, h1, s, cv, cp, w, ierr, herr, errormessagelength);
		//printf("T, D, P from PDFLSH   %10.4f,%10.4f,%10.4f\n", tt, d, p / 1000.0);

		////...General property calculation with inputs of p,h,x
		//PHFLSHdll(p, h, x, tt, dd, dl, dv, xliq, xvap, q, e, s, cv, cp, w, ierr, herr, errormessagelength);
		//printf("T, D, P from PHFLSH   %10.4f,%10.4f,%10.4f\n", tt, dd, p / 1000.0);

		////...General property calculation with inputs of p,s,x
		//PSFLSHdll(p, s, x, tt, dd, dl, dv, xliq, xvap, q, e, h1, cv, cp, w, ierr, herr, errormessagelength);
		//printf("T, D, P from PSFLSH   %10.4f,%10.4f,%10.4f\n", tt, dd, p / 1000.0);

		////...General property calculation with inputs of d,h,x
		//DHFLSHdll(d, h, x, tt, pp, dl, dv, xliq, xvap, q, e, s, cv, cp, w, ierr, herr, errormessagelength);
		//printf("T, D, P from DHFLSH   %10.4f,%10.4f,%10.4f\n", tt, d, pp / 1000.0);

		////...General property calculation with inputs of t,h,x
		////     kr--flag specifying desired root for multi-valued inputs:
		////         1=return lower density root
		////         2=return higher density root
		//long kr = 1;
		//THFLSHdll(t, h, x,
		//	kr, pp, dd, dl, dv, xliq, xvap, q, e, s, cv, cp, w, ierr, herr, errormessagelength);
		//printf("T, D, P from THFLSH   %10.4f,%10.4f,%10.4f\n", t, dd, pp / 1000.0);

		////...Other general property calculation routines are given in FLSH_SUB.FOR
		////...and FLASH2.FOR

		////...Calculate melting pressure
		//t = 100.0;
		//MELTTdll(t, x, p, ierr, herr, errormessagelength);
		//printf("Melting pressure(MPa) %10.4f,%10.4f\n", p / 1000.0, t);

		////...Calculate melting temperature
		//MELTPdll(p, x, tt, ierr, herr, errormessagelength);
		//printf("Melting temperature(K)%10.4f,%10.4f\n", tt, p / 1000.0);

		////...Calculate sublimation pressure
		//t = 200.0;
		//SUBLTdll(t, x, p, ierr, herr, errormessagelength);
		//printf("Sublimation pr.(kPa)  %10.4f,%10.4f\n", p, t);

		////...Calculate sublimation temperature
		//SUBLPdll(p, x, tt, ierr, herr, errormessagelength);
		//printf("Sublimation temp.(K)  %10.4f,%10.4f\n", tt, p);

		////...Get limits of the equations and check if t,d,p is a valid point
		////...Equation of state
		////     call LIMITK ('EOS',1,t,d,p,tmin,tmax,Dmax,pmax,ierr,herr)
		////...Viscosity equation
		////     call LIMITK ('ETA',1,t,d,p,tmin,tmax,Dmax,pmax,ierr,herr)
		////...Thermal conductivity equation
		////     call LIMITK ('TCX',1,t,d,p,tmin,tmax,Dmax,pmax,ierr,herr)

		////...Other routines are given in UTILITY.FOR

	}

	void setSubstance(char* fldname) 
	{

		//此处显式声明dll的位置，相对位置或者绝对位置都可以
		//Exlicitely set the fluid file PATH
		char *FLD_PATH;
		FLD_PATH = "fluids\\";
		//	  strcpy(hf,FLD_PATH);
		//	  strcpy(hfmix,FLD_PATH);

		/* ********************************************/
		//...initialize the program and set the pure fluid component name
		i = 1;
		strcpy(hf, fldname);
		// strcpy(hfmix,"hmx.bnc");
		strcpy(hrf, "DEF");
		strcpy(herr, "Ok");

		//...Call SETUP to initialize the program
		SETUPdll(i, hf, hfmix, hrf, ierr, herr,
			refpropcharlength*ncmax, refpropcharlength,
			lengthofreference, errormessagelength);
		if (ierr != 0) printf("%s\n", herr);

		double wm, ttp, tnbp, tc, pc, dc, zc, acf, dip, rgas;
		long info_index = 1;
		INFOdll(info_index, wm, ttp, tnbp, tc, pc, dc, zc, acf, dip, rgas);
		//printf("WM,ACF,DIP,TTP,TNBP   %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", wm, acf, dip, ttp, tnbp);
		//printf("TC,PC,DC,RGAS         %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", tc, pc, dc, rgas);

		printf("substance is selected.\n");

	}
	//根据T,P,计算物性的第一种方法,确定气态还是液态:
	//j = 1，液体；j = 2，气体
	//If phase is known : (j = 1: Liquid, j = 2 : Vapor)
	double* byTP(double m_T,double m_P,long j)
	{
		double t = m_T;
		double p = m_P;
		//printf("给定温度t = %f，压力 p = %f，x[0] = %f \n", t, p, x);
		//...Calculate d from t,p,x
		//...If phase is known: (j=1: Liquid, j=2: Vapor)
		double d, q, e, h, s, cv, cp, w, b, c,
			dpdrho, d2pdd2, dpdt, dhdt_d, dhdt_p, dhdp_t, dhdp_d,
			sigma, dhdd_t, dhdd_p, eta, tcx, pp, tt, hjt, h1, dd;
		long tmp_int = 0;
		TPRHOdll(t, p, x, j, tmp_int, d, ierr, herr, errormessagelength);
		//printf("T,P,D                 %10.4f,%10.4f,%10.4f\n", t, p, d);

		int size = 27;
		double* pdata = new double[size];
		
		pdata[0] = d;
		//pdata[1] = q;
		return pdata;

	}

	//根据T,P,计算物性的第二种方法, 不清楚气态还是液态
	double byTP2(double m_T, double m_P)
	{
		double t = m_T;
		double p = m_P;
		double dl, dv;
		double d, q, e, h, s, cv, cp, w, b, c,
			dpdrho, d2pdd2, dpdt, dhdt_d, dhdt_p, dhdp_t, dhdp_d,
			sigma, dhdd_t, dhdd_p, eta, tcx, pp, tt, hjt, h1, dd;
		TPFLSHdll(t, p, x, d, dl, dv, xliq, xvap, q, e, h, s, cv, cp, w, ierr, herr, errormessagelength);
		//printf("T,P,D,H,CP            %10.4f,%10.4f,%10.4f,%10.4f,%10.4f\n", t, p, d, h, cp);
		//printf("密度d= %f (kg/m3)", d);

		return d;

	}

private:

};

Substance::Substance()
{
	printf("constructor is called.\n");
	initialize();
}

Substance::~Substance()
{
}


int main()
{
	Substance Nitrogen;
	char* fldname = "nitrogen.fld";
	Nitrogen.setSubstance(fldname);
	Nitrogen.showSubstance();
	double* pd  = Nitrogen.byTP(77,2000,1);
	printf(" density is: %f(mol/L).\n",pd[0]);

	for (int px = 1000; px < 2000; px = px + 100)
	{
		pd = Nitrogen.byTP(77.0, px*1.0, 1);
		printf("density is 77 (K)，pressure is %d (KPa)， density = %f6.2 (mol/L).\n",px,pd[0]);
	}
	return 0;
}
