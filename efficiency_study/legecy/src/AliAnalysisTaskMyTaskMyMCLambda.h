#ifndef AliAnalysisTaskMyTaskMyMCLambda_cxx
#define AliAnalysisTaskMyTaskMyMCLambda_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliPIDResponse.h"
#include "AliESDEvent.h"
#include "AliEventCuts.h"

//class AliAnalysisUtils;

//#include "AliAnalysisTaskMy.h"

class AliAnalysisTaskMyTaskMyMCLambda : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskMyTaskMyMCLambda();
  AliAnalysisTaskMyTaskMyMCLambda(const char *name);
  virtual ~AliAnalysisTaskMyTaskMyMCLambda();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  bool IsGoodV0(AliAODv0 *aodV0);
  bool IsGoodDaughterTrack(const AliAODTrack *track);
  int  GetLambdaCode(const AliAODTrack *pTrack, const AliAODTrack *nTrack);
  double MyRapidity(double rE, double rPz) const;


 private:

  AliAODEvent *fAOD;    //! AOD object
  AliPIDResponse *fPID; //! PID object
  AliMultSelection *fMultSel; //! MultSelection object
  AliVVertex *fVtx;     //! vertexer object
  AliEventCuts fEventCuts; //! event cuts object
  AliMCEvent *fMC; //!

  double fVertex[3];
  TList* fOutputList;

  //hadron
  TH1D* h_eta_poshadron_mc;
  TH1D* h_eta_neghadron_mc;
  TH1D* h_eta_poshadron_rc;
  TH1D* h_eta_neghadron_rc;

  TH1D* h_pt_poshadron_mc;
  TH1D* h_pt_neghadron_mc;
  TH1D* h_pt_poshadron_rc;
  TH1D* h_pt_neghadron_rc;

  //pion
  TH1D* h_eta_pospion_mc;
  TH1D* h_eta_negpion_mc;
  TH1D* h_eta_pospion_rc;
  TH1D* h_eta_negpion_rc;
  TH1D* h_eta_pospion_rc_real;
  TH1D* h_eta_negpion_rc_real;


  TH1D* h_pt_pospion_mc;
  TH1D* h_pt_negpion_mc;
  TH1D* h_pt_pospion_rc;
  TH1D* h_pt_negpion_rc;
  TH1D* h_pt_pospion_rc_real;
  TH1D* h_pt_negpion_rc_real;

  //proton
  TH1D* h_eta_proton_mc;
  TH1D* h_eta_antiproton_mc;
  TH1D* h_eta_proton_rc;
  TH1D* h_eta_antiproton_rc;
  TH1D* h_eta_proton_rc_real;
  TH1D* h_eta_antiproton_rc_real;

  TH1D* h_pt_proton_mc;
  TH1D* h_pt_antiproton_mc;
  TH1D* h_pt_proton_rc;
  TH1D* h_pt_antiproton_rc;
  TH1D* h_pt_proton_rc_real;
  TH1D* h_pt_antiproton_rc_real;

  //lambda
  TH1D* h_rap_lambda_mc;
  TH1D* h_rap_antilambda_mc;
  TH1D* h_rap_lambda_rc;
  TH1D* h_rap_antilambda_rc;
  TH1D* h_rap_lambda_rc_real;
  TH1D* h_rap_antilambda_rc_real;

  TH1D* h_pt_lambda_mc;
  TH1D* h_pt_antilambda_mc;
  TH1D* h_pt_lambda_rc;
  TH1D* h_pt_antilambda_rc;
  TH1D* h_pt_lambda_rc_real;
  TH1D* h_pt_antilambda_rc_real;

  TH1D* h_mass_lambda_rc;
  TH1D* h_mass_antilambda_rc;
  TH1D* h_mass_lambda_rc_real;
  TH1D* h_mass_antilambda_rc_real;

  //proton dca
  TH2D* h_pt_dca_proton_origin;
  TH2D* h_pt_dca_proton_material;
  TH2D* h_pt_dca_proton_lambdadecay;
  TH2D* h_pt_dca_proton_otherdecay;

  //proton dca
  TH2D* h_pt_dca_antiproton_origin;
  TH2D* h_pt_dca_antiproton_material;
  TH2D* h_pt_dca_antiproton_lambdadecay;
  TH2D* h_pt_dca_antiproton_otherdecay;

  AliAnalysisTaskMyTaskMyMCLambda(const AliAnalysisTaskMyTaskMyMCLambda&);
  AliAnalysisTaskMyTaskMyMCLambda& operator=(const AliAnalysisTaskMyTaskMyMCLambda&);

  ClassDef(AliAnalysisTaskMyTaskMyMCLambda, 1);
};

#endif
