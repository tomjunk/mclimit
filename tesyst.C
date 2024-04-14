#include <iostream>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <stdio.h>
#include <stdlib.h>
#include <TCanvas.h>
#include <TFormula.h>
#include <TF1.h>
#include <TF2.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TMath.h>

// horizontal is the interpolation style which calls csm_pvmorph (or csm_pvmorph_2d)
// vertical interpolation just does a linear interpolation bin-by-bin, but not
// letting any bin go below zero

typedef enum {
  CSM_INTERP_HORIZONTAL,
  CSM_INTERP_VERTICAL,
  CSM_INTERP_HORIZONTAL_EXTRAP,
  CSM_INTERP_VERTICAL_EXTRAP
} INTERPSTYLE;


// interpolate histogram with a pvmorph-style procedure, or with vertical interpolation.
// csm_interpolate_histogram interpolates the errors too, while csm_interpolate_histogram_noerr
// is a speedup which interpolates only the bin contents and not the errors

void csm_interpolate_histogram(TH1*,Double_t,TH1*,Double_t,TH1*,Double_t,INTERPSTYLE);

void csm_interpolate_histogram_noerr(TH1*,Double_t,TH1*,Double_t,TH1*,Double_t,INTERPSTYLE);

// version to be used with cascading shape errors -- needs a central shape, a varied shape,
// and a shape to apply the variations to (which may not be either of the above, but the result
// of a different shape variation) startshape.  The output is outshape.

// csm_interpolate histogram2 calls csm_interpolate_histogram3 twice, once for the
// bin contents, once for the errors.

void
csm_interpolate_histogram2(TH1* central, Double_t paramcentral,
			   TH1* varied, Double_t paramvaried,
			   TH1* startshape, 
			   TH1* outshape,
			   Double_t param,
			   INTERPSTYLE istyle);

// here's a version that just interpolates the bin contents but not the errors
//  (twice as fast)

void
csm_interpolate_histogram2_noerr(TH1* central, Double_t paramcentral,
				 TH1* varied, Double_t paramvaried,
				 TH1* startshape, 
				 TH1* outshape,
				 Double_t param,
				 INTERPSTYLE istyle);

// this routine just interpolates the histogram contents

void
csm_interpolate_histogram3(TH1* central, Double_t paramcentral,
			   TH1* varied, Double_t paramvaried,
			   TH1* startshape, 
			   TH1* outshape,
			   Double_t param,
			   INTERPSTYLE istyle);

void csm_pvmc(Int_t nb, Double_t *dist1, Double_t *dist2, Double_t *dist3, Double_t *distn,
	      Double_t par1, Double_t par2, Double_t parn);

void csm_pvmc2d(Int_t nx, Int_t ny, Double_t *xydist1, 
                Double_t *xydist2, Double_t *xydist3, Double_t *xydistn,
                Double_t par1, Double_t par2, Double_t parn);

void csm_yproj(Int_t nx, Int_t ny, Double_t *xydist, Double_t *ydist);

void csm_ycont(Int_t ny, Double_t *ydist1, Double_t *ydist2,
               Double_t *ydist3, Double_t *ydistn,
               Double_t *alpha1, Double_t *alpha2, Double_t *alpha3);

void csm_ycontaux(Int_t ny, Double_t *y, Double_t *yn,
                  Double_t *alpha);

void csm_acnvec2(Double_t *vec, Int_t n);


#define PVMORPH_MAXBINS 5000


/*------------------------------------------------------------------------*/

// interpolate 1D histograms and 2D histograms
// histo a corresponds to parameter xa, histo b corresponds to xb.
// xc is input, and histogram c is the interpolated output

// new version -- rely on the more general interpolator with three inputs, but reduce the argument
// count for backward compatibility

// csm_interpolate_histogram interpolates the bin contents and errors

inline void
csm_interpolate_histogram(TH1* a, Double_t xa, 
			  TH1* b, Double_t xb,
			  TH1* c, Double_t xc,
			  INTERPSTYLE istyle)
{
  csm_interpolate_histogram2(a,xa,b,xb,a,c,xc,istyle);
}

// csm_interpolate_histogram_noerr interpolates just the bin contents but not the errors

void
csm_interpolate_histogram_noerr(TH1* a, Double_t xa, 
				TH1* b, Double_t xb,
				TH1* c, Double_t xc,
				INTERPSTYLE istyle)
{
  csm_interpolate_histogram2_noerr(a,xa,b,xb,a,c,xc,istyle);
}

// interpolate 1D histograms and 2D histograms 
// histo a corresponds to parameter xa, histo b corresponds to xb.
// xc is input, and histogram c is the interpolated output
// d is the histogram to apply the shift given by a and b to, for compounded interpolations.

// approximate attempt to interpolate the uncertainties too.  Problem is, an interpolated
// histogram is a long sum of pieces interpolated from the same central value histogram,
// and thus the errors are correlated in interesting ways.
// A subterfuge -- jut linearly interpolate the errors in the same way that the
// bin contents are linearly interpolated.  It's not a full error propagation.  Halfway interpolations
// really are averages of statistically uncertain histograms, and thus the error on the average should
// be a bit better than the error on either one.  But itnterpolate again, and correlations have to be
// taken into account to do it right.
// we've also lost at this point whether the errors need to be interpolated, but let's
// do them for all histograms.
// speedup 9 Dec 2007 -- avoid cloning TH1's as this is slow

inline void
csm_interpolate_histogram2(TH1* a, Double_t xa, 
			   TH1* b, Double_t xb,
			   TH1* d,
			   TH1* c, Double_t xc,
			   INTERPSTYLE istyle)
{
  Int_t i,j;
  Double_t xtmp;
  Int_t nbinsa = a->GetNbinsX();
  Int_t nbinsb = b->GetNbinsX();
  Int_t nbinsc = c->GetNbinsX();
  Int_t nbinsd = d->GetNbinsX();
  Int_t nbinsya = a->GetNbinsY();
  Int_t nbinsyb = b->GetNbinsY();
  Int_t nbinsyc = c->GetNbinsY();
  Int_t nbinsyd = d->GetNbinsY();

  if (nbinsa != nbinsb)
    {
      std::cout << "nbins mismatch1 in csm_interpolate_histogram2: " << nbinsa << " " << nbinsb << std::endl;
    }
  if (nbinsb != nbinsc)
    {
      std::cout << "nbins mismatch2 in csm_interpolate_histogram2: " << nbinsb << " " << nbinsc << std::endl;
    }
  if (nbinsc != nbinsd)
    {
      std::cout << "nbins mismatch3 in csm_interpolate_histogram2: " << nbinsc << " " << nbinsd << std::endl;
    }
  if (nbinsya != nbinsyb)
    {
      std::cout << "nbinsy mismatch1 in csm_interpolate_histogram2: " << nbinsya << " " << nbinsyb << std::endl;
    }
  if (nbinsyb != nbinsyc)
    {
      std::cout << "nbinsy mismatch2 in csm_interpolate_histogram2 " << nbinsyb << " " << nbinsyc << std::endl;
    }
  if (nbinsyc != nbinsyd)
    {
      std::cout << "nbinsy mismatch3 in csm_interpolate_histogram2: " << nbinsyc << " " << nbinsyd << std::endl;
    }

  if (xb == xa)
    {
      std::cout << "xb == xa in csm_interpolate_histogram2 " << xa << std::endl;
      std::cout << "fatal error -- exiting." << std::endl;
      exit(0);
    }

  // interpolate contents

  csm_interpolate_histogram3(a,xa,b,xb,d,c,xc,istyle);

  // swap errors and contents and interpolate again  (approximate method for evaluating
  // errors on interpolated histograms)
  // be careful to swap only once, even if some pointers are repeated.

  if (nbinsya == 1)
    {
      for (i=1;i<=nbinsa;i++)
	{
	  xtmp = a->GetBinContent(i);
	  a->SetBinContent(i,a->GetBinError(i));
	  a->SetBinError(i,xtmp);
	  if (a != b)
	    {
	      xtmp = b->GetBinContent(i);
	      b->SetBinContent(i,b->GetBinError(i));
	      b->SetBinError(i,xtmp);
	    }
	  // c is the output histogram -- hopefully it is not the same as one of the input histograms
          xtmp = c->GetBinContent(i);
          // c->SetBinContent(i,c->GetBinError(i));
	  c->SetBinError(i,xtmp);

	  if (a != d && b != d)
	    {
	      xtmp = d->GetBinContent(i);
	      d->SetBinContent(i,d->GetBinError(i));
	      d->SetBinError(i,xtmp);
	    }
	}
    }
  else
    {
      for (i=1;i<=nbinsa;i++)
	{
	  for (j=1;j<=nbinsya;j++)
	    {
	      xtmp = a->GetBinContent(i,j);
	      a->SetBinContent(i,j,a->GetBinError(i,j));
	      a->SetBinError(i,j,xtmp);
	      if (a != b)
		{
		  xtmp = b->GetBinContent(i,j);
		  b->SetBinContent(i,j,b->GetBinError(i,j));
		  b->SetBinError(i,j,xtmp);
		}
	      xtmp = c->GetBinContent(i,j);
	      //c->SetBinContent(i,j,c->GetBinError(i,j));
	      c->SetBinError(i,j,xtmp);
	      if (a != d && b != d)
		{
		  xtmp = d->GetBinContent(i,j);
		  d->SetBinContent(i,j,d->GetBinError(i,j));
		  d->SetBinError(i,j,xtmp);
		}
	    }
	}
    }

  // interpolate the errors now and swap them back -- put the
  // original histograms back together again too

  csm_interpolate_histogram3(a,xa,b,xb,d,c,xc,istyle);

  if (nbinsya == 1)
    {
      for (i=1;i<=nbinsa;i++)
	{
	  xtmp = a->GetBinContent(i);
	  a->SetBinContent(i,a->GetBinError(i));
	  a->SetBinError(i,xtmp);
	  if (a != b)
	    {
	      xtmp = b->GetBinContent(i);
	      b->SetBinContent(i,b->GetBinError(i));
	      b->SetBinError(i,xtmp);
	    }
	  xtmp = c->GetBinContent(i);
	  c->SetBinContent(i,c->GetBinError(i));
	  c->SetBinError(i,xtmp);
	  if (a != d && b != d)
	    {
	      xtmp = d->GetBinContent(i);
	      d->SetBinContent(i,d->GetBinError(i));
	      d->SetBinError(i,xtmp);
	    }
	}
    }
  else
    {
      for (i=1;i<=nbinsa;i++)
	{
	  for (j=1;j<=nbinsya;j++)
	    {
	      xtmp = a->GetBinContent(i,j);
	      a->SetBinContent(i,j,a->GetBinError(i,j));
	      a->SetBinError(i,j,xtmp);
	      if (a != b)
		{
		  xtmp = b->GetBinContent(i,j);
		  b->SetBinContent(i,j,b->GetBinError(i,j));
		  b->SetBinError(i,j,xtmp);
		}
	      xtmp = c->GetBinContent(i,j);
	      c->SetBinContent(i,j,c->GetBinError(i,j));
	      c->SetBinError(i,j,xtmp);
	      if (a != d && b != d)
		{
		  xtmp = d->GetBinContent(i,j);
		  d->SetBinContent(i,j,d->GetBinError(i,j));
		  d->SetBinError(i,j,xtmp);
		}
	    }
	}
    }
}

void
csm_interpolate_histogram2_noerr(TH1* a, Double_t xa, 
				 TH1* b, Double_t xb,
				 TH1* d,
				 TH1* c, Double_t xc,
				 INTERPSTYLE istyle)
{
  Int_t nbinsa = a->GetNbinsX();
  Int_t nbinsb = b->GetNbinsX();
  Int_t nbinsc = c->GetNbinsX();
  Int_t nbinsd = d->GetNbinsX();
  Int_t nbinsya = a->GetNbinsY();
  Int_t nbinsyb = b->GetNbinsY();
  Int_t nbinsyc = c->GetNbinsY();
  Int_t nbinsyd = d->GetNbinsY();

  if (nbinsa != nbinsb)
    {
      std::cout << "nbins mismatch1 in csm_interpolate_histogram2_noerr: " << nbinsa << " " << nbinsb << std::endl;
    }
  if (nbinsb != nbinsc)
    {
      std::cout << "nbins mismatch2 in csm_interpolate_histogram2_noerr: " << nbinsb << " " << nbinsc << std::endl;
    }
  if (nbinsc != nbinsd)
    {
      std::cout << "nbins mismatch3 in csm_interpolate_histogram2_noerr: " << nbinsc << " " << nbinsd << std::endl;
    }
  if (nbinsya != nbinsyb)
    {
      std::cout << "nbinsy mismatch1 in csm_interpolate_histogram2_noerr: " << nbinsya << " " << nbinsyb << std::endl;
    }
  if (nbinsyb != nbinsyc)
    {
      std::cout << "nbinsy mismatch2 in csm_interpolate_histogram2_noerr " << nbinsyb << " " << nbinsyc << std::endl;
    }
  if (nbinsyc != nbinsyd)
    {
      std::cout << "nbinsy mismatch3 in csm_interpolate_histogram2_noerr: " << nbinsyc << " " << nbinsyd << std::endl;
    }

  if (xb == xa)
    {
      std::cout << "xb == xa in csm_interpolate_histogram2_noerr " << xa << std::endl;
      std::cout << "fatal error -- exiting." << std::endl;
      exit(0);
    }

  // interpolate just the bin contents

  csm_interpolate_histogram3(a,xa,b,xb,d,c,xc,istyle);

}


inline void
csm_interpolate_histogram3(TH1* a, Double_t xa, 
			   TH1* b, Double_t xb,
			   TH1* d,
			   TH1* c, Double_t xc,
			   INTERPSTYLE istyle)
{
  Double_t hnorma,hnormb,hnormc,hnormd,hnormci;
  Int_t i,j;
  Double_t gbc;

  Int_t nbinsa = a->GetNbinsX();
  Int_t nbinsb = b->GetNbinsX();
  Int_t nbinsc = c->GetNbinsX();
  Int_t nbinsd = d->GetNbinsX();
  Int_t nbinsya = a->GetNbinsY();
  Int_t nbinsyb = b->GetNbinsY();
  Int_t nbinsyc = c->GetNbinsY();
  Int_t nbinsyd = d->GetNbinsY();

  if (a->Integral()<=0 || b->Integral()<=0)
    { 
      for (i=1;i<=nbinsc;i++)
	{
	  for (j=1;j<=nbinsyc;j++)
	    { c->SetBinContent(i,j,0);
	    }
	} 
      //c->Reset();
      return;
    }
    
  if (nbinsya == 1)
    {
      Double_t *dista = new Double_t[nbinsa];
      Double_t *distb = new Double_t[nbinsb];
      Double_t *distc = new Double_t[nbinsc];
      Double_t *distd = new Double_t[nbinsd];

      hnorma = 0;
      hnormb = 0;
      hnormd = 0;
      for (i=0;i<nbinsa;i++)
        { dista[i] = a->GetBinContent(i+1); hnorma += dista[i]; }
      for (i=0;i<nbinsb;i++)
        { distb[i] = b->GetBinContent(i+1); hnormb += distb[i]; }
      for (i=0;i<nbinsd;i++)
        { distd[i] = d->GetBinContent(i+1); hnormd += distd[i]; }

      hnormc = hnorma + (xc-xa)*(hnormb-hnorma)/(xb-xa);
      // linearly interpolate the normalization between the central value and
      // the varied template.
      hnormc = hnormd*(hnormc/hnorma); // scale the normalization with the new template
      if (hnormc<0) hnormc = 0;

      if (istyle == CSM_INTERP_HORIZONTAL || istyle == CSM_INTERP_HORIZONTAL_EXTRAP)
	{
	  csm_pvmc(nbinsa,dista,distb,distd,distc,xa,xb,xc);
	  hnormci = 0;
	  for (i=0;i<nbinsc;i++) { hnormci += distc[i]; }

	  for (i=0;i<nbinsc;i++)
	    {
	      if (hnormci != 0)
		{
		  c->SetBinContent(i+1,distc[i]*hnormc/hnormci);
		}
	      else 
		{
		  c->SetBinContent(i+1,0);
		}
	    }
	}
      else if (istyle == CSM_INTERP_VERTICAL || istyle == CSM_INTERP_VERTICAL_EXTRAP)
	{
	  for (i=0;i<nbinsa;i++)
	    {
	      gbc = distd[i] + ((xc-xa)/(xb-xa))*(distb[i]-dista[i]);
	      if (gbc < 0) 
		{
		  gbc = 0;
		}

	      c->SetBinContent(i+1,gbc);
	    }
	}
      else
	{
	  std::cout << "csm_interpolate_histogram: unknown interpolation style " << istyle << std::endl;
	  exit(0);
	}

      //std::cout << xa << " " << xb << " " << xc << std::endl;

      delete[] dista;
      delete[] distb;
      delete[] distc;
      delete[] distd;
    }
  else         // 2d case
    {
      Double_t *distxya = new Double_t[nbinsa*nbinsya];
      Double_t *distxyb = new Double_t[nbinsb*nbinsyb];
      Double_t *distxyc = new Double_t[nbinsc*nbinsyc];
      Double_t *distxyd = new Double_t[nbinsd*nbinsyd];

      hnorma = 0;
      for (j=0;j<nbinsya;j++)
	{
	  for (i=0;i<nbinsa;i++)
	    {
	      gbc = a->GetBinContent(i+1,j+1);
	      distxya[i+nbinsa*j] = gbc;
	      hnorma += gbc;
	    }
	}

      hnormb = 0;
      for (j=0;j<nbinsyb;j++)
	{
	  for (i=0;i<nbinsb;i++)
	    {
	      gbc = b->GetBinContent(i+1,j+1);
	      distxyb[i+nbinsb*j] = gbc;
	      hnormb += gbc;
	    }
	}

      hnormd = 0;
      for (j=0;j<nbinsyd;j++)
	{
	  for (i=0;i<nbinsd;i++)
	    {
	      gbc = d->GetBinContent(i+1,j+1);
	      distxyd[i+nbinsb*j] = gbc;
	      hnormd += gbc;
	    }
	}

      hnormc = hnorma + (xc-xa)*(hnormb-hnorma)/(xb-xa);
      // linearly interpolate the normalization between the central value and
      // the varied template.
      hnormc = hnormd*(hnormc/hnorma); // scale the normalization with the new template
      if (hnormc<0) hnormc = 0;

      if (istyle == CSM_INTERP_HORIZONTAL || istyle == CSM_INTERP_HORIZONTAL_EXTRAP)
	{
          csm_pvmc2d(nbinsa,nbinsya,
                     distxya,
                     distxyb,
                     distxyd,
		     distxyc,
                     xa, xb, xc);

          hnormci = 0;
          for (j=0;j<nbinsyc;j++)
	    {
              for (i=0;i<nbinsc;i++)
	        {
	          hnormci += distxyc[i+nbinsc*j];
	        }
	    }
          for (j=0;j<nbinsyc;j++)
	    {
              for (i=0;i<nbinsc;i++)
	        {
		  if (hnormci != 0)
		    {
	              c->SetBinContent(i+1,j+1,distxyc[i+nbinsc*j]*hnormc/hnormci);
		    }
		  else
		    {
	              c->SetBinContent(i+1,j+1,0);
		    }
	        }
	    }
	}
      else if (istyle == CSM_INTERP_VERTICAL || istyle == CSM_INTERP_VERTICAL_EXTRAP)
	{
          for (j=0;j<nbinsyc;j++)
	    {
              for (i=0;i<nbinsc;i++)
	        {
		  gbc = distxyd[i+nbinsc*j] + ((xc-xa)/(xb-xa))*(distxyb[i+nbinsc*j]-distxya[i+nbinsc*j]);
		  if (gbc < 0)
		    {
		      gbc = 0;
		    }
	          c->SetBinContent(i+1,j+1,gbc);
	        }
	    }
	}
      else
	{
	  std::cout << "csm_interpolate_histogram: unknown interpolation style " << istyle << std::endl;
	  exit(0);
	}

      delete[] distxya;
      delete[] distxyb;
      delete[] distxyc;
      delete[] distxyd;
    }

}

/*------------------------------------------------------------------------*/

/* compounded interpolation -- dist1 = central value shape, dist2 = syst. varied shape,
   dist3 = shape to distort (may be the result of previous distortions for compounded shape
   variations), distn = resultant shape.  par1 = value of parameter for dist1.  par2 = value of
   parameter (like # of sigma) for dist2.  parn = value of parameter for the output histogram
   Built on the idea of d_pvmorph, but generalized a bit.  Returns a null histogram if any of
   the three input histograms has zero or negative total sum.
*/

//#define DEBUGPVMC

void csm_pvmc(Int_t nb, Double_t *dist1, Double_t *dist2, Double_t *dist3, Double_t *distn,
	      Double_t par1, Double_t par2, Double_t parn)
{
  Int_t nb3 = nb*3 + 3;
  Int_t i,j,k,k1,k2;
  Double_t total;
  Double_t wta,wtb;
  Double_t yd[nb3];
  Double_t id[nb3];
  Double_t xd[nb3];
  Double_t xdis[nb3];
  Double_t ydis[nb3];
  Double_t ydi[nb + 1];
  Int_t idx[nb3];
  Int_t ifirst;
  Double_t x1l,x2l,x3l,y1l,y2l,y3l;
  Double_t x1,x2,x3,y1,y2,y3;
  Double_t xloc,yloc,x1i,x2i,x3i;

  // default output -- empty distribution

  for (i=0;i<nb;i++)
    {
      distn[i] = 0.0;
    }

  // default index list

  for (i=0;i<nb3;i++)
    { idx[i] = i;}

  // parameter weights

  if (par2 != par1) 
    {
      wta = 1. - (parn-par1)/(par2-par1);
      wtb = 1. + (parn-par2)/(par2-par1);
    }
  else
    {
      wta = 0.5;
      wtb = 0.5;
    }

  // suppress warning messages in case of extrapolations

  //  if ( (parn>par1 && parn>par2) || (parn<par1 && parn<par2) )
  //  {
  //    cout << "CSM_PVMC: Histogram Extrapolation: " << parn << 
  //            " is not between " << par1 << " and " << par2 << endl;
  //  }

  // Fill cumulative distribution arrays -- squeeze out repeated entries
  // due to empty bins

  // The first point in the cumulative distributions has zero integral and
  // starts at the left-hand edge of the first bin with any value in it.
  // The id array says which distribution it came from, and the
  // xd array gives the x value at which the cumulative distribution is evaluated
  // (at the right-hand edge of the bin)

  j = 0;
  total = 0;
  for (i=0;i<nb;i++)
    {
      if (dist1[i] < 0)
	{ 
	  cout << "Negative bin entry found in dist1 in csm_pvmc" << endl;
	  cout << i << " " << dist1[i] << endl;
	  exit(0);
	}
      total += dist1[i];
    }
  if (total <= 0) return;

  yd[j] = 0;
  id[j] = 1;
  j++;
  ifirst = 1;
  for (i=0;i<nb;i++)
    {
      if (dist1[i] > 0)
	{ 
	  if (ifirst==1)
	    {
	      ifirst = 0;
	      xd[j-1] = (Double_t) i;
	    }
	  yd[j] = yd[j-1] + dist1[i]/total;
	  id[j] = 1;
	  xd[j] = (Double_t) i+1;
	  j++;
	}
    }

  total = 0;
  for (i=0;i<nb;i++)
    {
      if (dist2[i] < 0)
	{ 
	  cout << "Negative bin entry found in dist2 in csm_pvmc" << endl;
	  cout << i << " " << dist2[i] << endl;
	  exit(0);
	}
      total += dist2[i];
    }
  if (total <= 0) return;
  yd[j] = 0;
  id[j] = 2;
  j++;
  ifirst = 1;
  for (i=0;i<nb;i++)
    {
      if (dist2[i]>0)
	{
	  if (ifirst==1)
	    {
	      ifirst = 0;
	      xd[j-1] = (Double_t) i;
	    }
          yd[j] = yd[j-1] + dist2[i]/total;
	  id[j] = 2;
	  xd[j] = (Double_t) i+1;
	  j++;
	}
    }

  total = 0;
  for (i=0;i<nb;i++)
    {
      if (dist3[i] < 0)
	{ 
	  cout << "Negative bin entry found in dist3 in csm_pvmc" << endl;
	  cout << i << " " << dist3[i] << endl;
	  exit(0);
	}
      total += dist3[i];
    }
  if (total <= 0) return;
  yd[j] = 0;
  id[j] = 3;
  j++;
  ifirst = 1;
  for (i=0;i<nb;i++)
    {
      if (dist3[i]>0)
	{
	  if (ifirst==1)
	    {
	      ifirst = 0;
	      xd[j-1] = (Double_t) i;
	    }
          yd[j] = yd[j-1] + dist3[i]/total;
	  id[j] = 3;
	  xd[j] = (Double_t) i+1;
	  j++;
	}
    }

  // Sort all of the edges of the cumulative distribution functions
  // j is the number of entries in the yd, xd and id arrays

  TMath::Sort(j,yd,idx,0);

#ifdef DEBUGPVMC
  for (i=0;i<j;i++)
    {
      cout << i << " " << xd[i] << " " << " " << yd[i] << " " << id[i] << endl;
    }
  cout << "Sort index" << endl;
  for (i=0;i<j;i++)
    {
      cout << i << " " << idx[i] << endl;
    }
#endif

  x1l = 0;
  x2l = 0;
  x3l = 0;

  y1l = 0;
  y2l = 0;
  y3l = 0;

  x1 = 0;
  x2 = 0;
  x3 = 0;

  y1 = 0;
  y2 = 0;
  y3 = 0;

  // the three lowest points in the sort should all have zero integral --
  // interpolate the x's of these

  for (i=0;i<3;i++)
    {
      if ( id[idx[i]] == 1 )
	{
	  x1 = xd[idx[i]];
	  y1 = yd[idx[i]]; // should be zero
	}
      else if ( id[idx[i]] == 2 )
	{
	  x2 = xd[idx[i]];
	  y2 = yd[idx[i]];  // should be zero
	}
      else if ( id[idx[i]] == 3 )
	{
	  x3 = xd[idx[i]];
	  y3 = yd[idx[i]];  // should be zero
	}
    }
  // don't have the other ends of the line segments yet -- find them as we go along.

#ifdef DEBUGPVMC
  cout << "first bins: " << x1 << " " << x2 << " " << x3 << endl;
#endif
  y1l = y1;
  y2l = y2;
  y3l = y3;
  x1l = x1;
  x2l = x2;
  x3l = x3;

  // first point on interpolated curve -- zero integral.

  k = 0;
  xdis[k] = wta*x1l + wtb*x2l - x1l + x3l;
  xdis[k] = min((Double_t) (nb+1),max(0.0,xdis[k]));
  ydis[k] = 0;

#ifdef DEBUGPVMC
  cout << "first point: " << xdis[0] << " " << ydis[0] << endl;
#endif

  for (i=3;i<j;i++)
    {
      xloc = xd[idx[i]];
      yloc = yd[idx[i]];

      if (id[idx[i]] == 1 )
	{
	  x1l = x1;
	  y1l = y1;
	  x1 = xloc;
	  y1 = yloc;

	  if (yloc>y2)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==2)
		    {
		      y2l = y2;
		      x2l = x2;
		      y2 = yd[idx[k1]];
		      x2 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	  if (yloc>y3)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==3)
		    {
		      y3l = y3;
		      x3l = x3;
		      y3 = yd[idx[k1]];
		      x3 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	}
      else if (id[idx[i]] == 2 )
	{
	  x2l = x2;
	  y2l = y2;
	  x2 = xloc;
	  y2 = yloc;

	  if (yloc>y1)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==1)
		    {
		      y1l = y1;
		      x1l = x1;
		      y1 = yd[idx[k1]];
		      x1 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	  if (yloc>y3)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==3)
		    {
		      y3l = y3;
		      x3l = x3;
		      y3 = yd[idx[k1]];
		      x3 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	}
      else if (id[idx[i]] == 3 )
	{
	  x3l = x3;
	  y3l = y3;
	  x3 = xloc;
	  y3 = yloc;

	  if (yloc>y2)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==2)
		    {
		      y2l = y2;
		      x2l = x2;
		      y2 = yd[idx[k1]];
		      x2 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	  if (yloc>y1)
	    {
	      for (k1=i+1;k1<j;k1++)
		{
		  if (id[idx[k1]]==1)
		    {
		      y1l = y1;
		      x1l = x1;
		      y1 = yd[idx[k1]];
		      x1 = xd[idx[k1]];
		      break;
		    }
		}
	    }
	}

      if (yloc>ydis[k] && ydis[k] < 0.999999999)
        {

#ifdef DEBUGPVMC
	  cout << "Interpolating: " << x1 << " " << x2 << " " << x3 << endl;
	  cout << "Interpolating: " << x1l << " " << x2l << " " << x3l << endl;
	  cout << "Interpolating: " << y1 << " " << y2 << " " << y3 << endl;
	  cout << "Interpolating: " << y1l << " " << y2l << " " << y3l << endl;
#endif

	  k++;
	  ydis[k] = yloc;
	  if (yloc == y1l)
	    {
	      x1i = x1l;
	    }
	  else
	    {
	      x1i = x1l + (yloc-y1l)*(x1-x1l)/(y1-y1l);
	    }
	  if (yloc == y2l)
	    {
	      x2i = x2l;
	    }
	  else
	    {
	      x2i = x2l + (yloc-y2l)*(x2-x2l)/(y2-y2l);
	    }
	  if (yloc == y3l)
	    {
	      x3i = x3l;
	    }
	  else
	    {
	      x3i = x3l + (yloc-y3l)*(x3-x3l)/(y3-y3l);
	    }
	  xdis[k] = x3i + wta*x1i + wtb*x2i - x1i;
          xdis[k] = min((Double_t) (nb+1),max(0.0,xdis[k]));
	  if (xdis[k]<xdis[k-1])
	    {
	      k--;
	      ydis[k] = yloc;
	    }

#ifdef DEBUGPVMC
	  cout << "point: " << k << endl;
	  cout << "x1i, x2i, x3i: " << x1i << " " << x2i << " " << x3i << endl;
	  cout << "Interpolated: " << xdis[k] << " " << yloc << endl;
#endif
	}
    }

#ifdef DEBUGPVMC
  for (i=0;i<=k;i++)
    {
      cout << "IC before bin: " << i << " " << xdis[i] << " " << ydis[i] << endl;
    }
#endif


  // k is the index of the last entry in the xdis, ydis interpolated array.
  // find the places where the piecewise linear interpolated cumulative distribution
  // crosses the bin edges  the index on ydi is the low bin edge.

  ydi[0] = 0.0;
  for (i=1;i<(1 + nb);i++)
    {
      ydi[i] = 1.0;
    }

  Int_t k2last;
  k1 = 0;
  for (i=0;i<=k;i++)
    {
      k2last = k1;
      for (k2=k1+1;k2<(nb+1);k2++)
	{
	  if ( (Double_t) k2 < xdis[i] )
	    {
	      if (i==0)
		{
		  ydi[k2] = 0;
		}
	      else
		{
		  ydi[k2] = ydis[i-1] + ( (Double_t) k2  - xdis[i-1] )*
		    (ydis[i]-ydis[i-1])/(xdis[i]-xdis[i-1]);
#ifdef DEBUGPVMC
		  cout << "filling bins: " << i << " " << k1 << " " << k2 << " " << ydi[k2] << endl; 
#endif

		}
	      k2last = k2;
	    } 
	  if ( (Double_t) k2 > xdis[i] ) break;
	}
      k1 = k2last;
    }

#ifdef DEBUGPVMC
  for (i=0;i<(nb+1);i++)
    {
      cout << "interp. cumulative: " << i << " " << ydi[i] << endl;
    }
#endif

  // differentiate to get the output distn

  for (i=0;i<(nb);i++)
    {
      distn[i] = ydi[i+1] - ydi[i];
      if (distn[i]<0) 
	{
	  distn[i] = 0;
	  //cout << "negative element in interpolated histo: " << i << " " << distn[i] << endl;
	  //exit(0);
	}
    }
}

/*------------------------------------------------------------------------*/


/*  
    Re-coded version of d_pvmorph_2d from Alex Read.  C version from Tom Junk
    Added feature of compounding shape variations as systematic errors.
    February 2007

    ......Do a linear interpolation between three two-dimensional
    probability distributions (scatterplots) as a function
    of the characteristic parameter of the distribution, for use
    in both interpolation and in application of systematic uncertainties.
    xydist1 is the "central value" histogram
    xydist2 is the "systematically varied" histogram
    xydist3 is the histogram to apply the variation to
    xydistn is the output histogram.  See csm_pvmc
    for compounded systematic variation applicaiton in 1D
    (used repeatedly in here).

    This is a generalization of csm_pvmc. The 2d distribution
    can be move around and be stretched or squeezed in two
    dimenions but finite rotations (changes in the correlation)
    are poorly approximated.

    nx        : Number of x-bins in the input and output distributions.
    ny        : Number of y-bins in the input and output distributions.
    xydist1,xydist2,xydist3,xydistn
    : Bin contents of scatterplots. The arrays should be
    packed with the index running fastest over the x
    dimension.
    Contents are in xydist[ix+nx*iy]
    par1,par2,parn     : Values of the linear parameters that characterise the
    histograms (e.g. the Higgs mass).

    Output: xydistn.  Same binning as xydist1,xydist2,xydist3.
    Its memory must be allocated outside of the routine

*/

//#define DEBUGPVMC2D

void csm_pvmc2d(Int_t nx, Int_t ny, Double_t *xydist1, 
                Double_t *xydist2, Double_t *xydist3, Double_t *xydistn,
                Double_t par1, Double_t par2, Double_t parn)
{
  Double_t ydist1[ny],ydist2[ny],ydist3[ny],ydistn[ny];
  Double_t xtemp1[nx],xtemp2[nx],xtemp3[nx],xtempn[nx];
  Double_t alpha1[ny*ny],alpha2[ny*ny],alpha3[ny*ny];
  Int_t i,j,k;

  // Project xydist1,2,3 onto the y axis and normalize

  csm_yproj(nx,ny,xydist1,ydist1);
  csm_yproj(nx,ny,xydist2,ydist2);
  csm_yproj(nx,ny,xydist3,ydist3);

  // Interpolate the y-projections 

  csm_pvmc(ny,ydist1,ydist2,ydist3,ydistn,par1,par2,parn);

#ifdef DEBUGPVMC2D
  for (i=0;i<ny;i++)
    {
      cout << "iy: " << i << " " << ydist1[i] << " " << ydist2[i] << " " << ydistn[i] << endl;
    }
#endif

  // Find out which y bins of histograms 1,2,3 contribute
  // to each bin of the interpolated distribution ydistn

  csm_ycont(ny,ydist1,ydist2,ydist3,ydistn,alpha1,alpha2,alpha3);

  // Extract the x-distributions in the y-slice determined above
  // and interpolate them

  for (i=0;i<ny;i++)  // loop over resulting bins
    {
      for (k=0;k<nx;k++) xtemp1[k] = 0;
      for (k=0;k<nx;k++) xtemp2[k] = 0;
      for (k=0;k<nx;k++) xtemp3[k] = 0;
      for (j=0;j<ny;j++) // loop over contributing bins
	{
	  for (k=0;k<nx;k++)
	    {
	      xtemp1[k] += alpha1[j+ny*i]*xydist1[k+nx*j];
	      xtemp2[k] += alpha2[j+ny*i]*xydist2[k+nx*j];
	      xtemp3[k] += alpha3[j+ny*i]*xydist3[k+nx*j];
	    }
	}
      // Interpolate the x distributions
      csm_pvmc(nx,xtemp1,xtemp2,xtemp3,xtempn,par1,par2,parn);

      // Insert the interpolated x distribution into the final output dist
      for (k=0;k<nx;k++) xydistn[k+nx*i] = xtempn[k]*ydistn[i];
    }
}

/*
  Re-coded version of d_ypvscat from Alex Read.  C version from Tom Junk
  Project a scatterplot onto the y-axis.The
  projection is normalized so that the sum of the bin contents is 1.0.

  nx,ny    : Number of bins in the scatterplot for the x and y coordindates.
  The projection is done onto <ny> bins.
  xydist   : The 2-dimensional array of the probabilities
  ydist    : The 1-dimensional array of the 2d probabilities projected onto
  the y-axis.

  Inputs : nx,ny,xydist
  Outputs: ydist (ny is unchanged from input to output)
*/

void csm_yproj(Int_t nx, Int_t ny, Double_t *xydist, Double_t *ydist)
{
  Int_t i,j;
  Double_t total;

  for (i=0;i<ny;i++) ydist[i] = 0;
  total = 0;

  for (i=0;i<ny;i++)
    {
      for (j=0;j<nx;j++)
	{
          ydist[i] += xydist[j+nx*i];
	}
      total += ydist[i]; 
    }

  if (total>0)
    {
      for (i=0;i<ny;i++) ydist[i] /= total;
    }
}

/*
  Recoded d_getycont -- original by Alex Read, recoded by Tom Junk
  February 2007

  <ydist1> and <ydist2> and <ydist3>
  are the projections on the y-axis of three
  scatterplots which are going to be interpolated. <ydistn> is
  the interpolated 1d distribution which represent the projection
  of the interpolated scatterplot on the y-axis. This routine determines
  which bins of <ydist1> and <ydist2> and <ydist3>
  contribute and by what amount to
  each bin of <ydistn>. This information is used in csm_pvmc2d to 
  determine the input distributions in the x-direction of each
  y-bin: these are then interpolated and accumulated in the interpolated
  2d distribution.

  Inputs : ny,ydist1,ydist2,ydist3,ydistn
  Outputs: alpha1,alpha2,alpha3

  alpha1[iyc+ny*iy] encodes the contribution of bin iyc in ydist1
  to to bin iy in ydistn

*/

void csm_ycont(Int_t ny, Double_t *ydist1, Double_t *ydist2,
               Double_t *ydist3, Double_t *ydistn,
               Double_t *alpha1, Double_t *alpha2, Double_t *alpha3)
{
  Double_t y[ny+1];
  Double_t yn[ny+1];
  Int_t i;

  // Make arrays to describe the straight-line approximations
  // to the four cumulative distributions, y1,y2,y3,yn 
  // Make sure to start out with a point at 0

  yn[0] = 0;
  for (i=0;i<ny;i++) yn[i+1] = ydistn[i];
  csm_acnvec2(yn,ny+1);

  y[0] = 0;
  for (i=0;i<ny;i++) y[i+1] = ydist1[i];
  csm_acnvec2(y,ny+1);
  csm_ycontaux(ny,y,yn,alpha1);
#ifdef DEBUGPVMC2D
  for (i=0;i<ny+1;i++)
    {
      cout << "getting alpha1: " << i << " " << y[i] << " " << yn[i] << endl;
    }
  Int_t j; 
  for (i=0;i<ny;i++)
    {
      for (j=0;j<ny;j++)
	{
	  cout << i << " " << j << " " << alpha1[i+ny*j] << endl;
	}
    }
  
#endif

  y[0] = 0;
  for (i=0;i<ny;i++) y[i+1] = ydist2[i];
  csm_acnvec2(y,ny+1);
  csm_ycontaux(ny,y,yn,alpha2);

  y[0] = 0;
  for (i=0;i<ny;i++) y[i+1] = ydist3[i];
  csm_acnvec2(y,ny+1);
  csm_ycontaux(ny,y,yn,alpha3);
}

void csm_ycontaux(Int_t ny, Double_t *y, Double_t *yn,
                  Double_t *alpha)
{
  Int_t ny2;
  Int_t i,j;

  ny2 = ny*ny;

  // clear out the alpha array

  for (i=0;i<ny2;i++) alpha[i] = 0;

  // loop over bins and see what fraction each contributes

  for (i=0;i<ny;i++) // interpolated histogram
    {
      for (j=0;j<ny;j++) // contributing histogram bin
	{
          if (y[j+1]-y[j]>0)
	    {
	      // first case -- contributing bin entirely contained
	      // within the interpolated output bin.
	      if (y[j]>=yn[i] && y[j]<yn[i+1] && 
		  y[j+1]>=yn[i] && y[j+1]<yn[i+1])
		{
		  alpha[j+ny*i] = 1;
		}
	      // second case -- interpolated output bin is entirely
	      // contained within the contributing bin
	      else if (y[j]<yn[i] && y[j+1] >= yn[i+1])
		{
		  alpha[j+ny*i] = (yn[i+1]-yn[i])/(y[j+1]-y[j]);
		}
	      // third case -- contributing bin straddles the
	      // left edge of the interpolated output bin but ends inside
	      // the bin
	      else if (y[j]<yn[i] && y[j+1]>=yn[i] && y[j+1]<yn[i+1])
		{
		  alpha[j+ny*i] = (y[j+1]-yn[i])/(y[j+1]-y[j]);
		}
	      // fourth case -- contributing bin straddles the
	      // right edge of the interpolated output bin but starts inside
	      // the output bin
	      else if (y[j]>=yn[i] && y[j]<yn[i+1] && y[j+1]>=yn[i+1])
		{
		  alpha[j+ny*i] = (yn[i+1]-y[j])/(y[j+1]-y[j]);
		}
	      // non-overlapping case -- do nothing.
	      // save some time if we're beyond the edge
	      if (y[j]>yn[i+1]) break; 
	    }
	}
    } 
}

// Integrate and normalize a vector -- pretty much the same as
// csm_acnvec

void csm_acnvec2(Double_t *vec, Int_t n)
{
  Int_t i;
  Double_t tot;
  for (i=1;i<n;i++)
    {
      vec[i] += vec[i-1];
    }
  tot = vec[n-1];
  if (tot > 0)
    {
      for (i=0;i<n;i++) vec[i] /= tot;
    }
}

void te(int iswitch=0)
{
  TH1F hleft("hleft"," ",100,0,1);
  TH1F hright("hright"," ",100,0,1);
  TH1F hinterph("hinterph"," ",100,0,1);
  TH1F hinterpv("hinterpv"," ",100,0,1);

  TF1 form1("form1","exp(-(x-0.2)*(x-0.2)/0.005)");
  TF1 form2("form2","exp(-(x-0.7)*(x-0.7)/0.01)");

  if (iswitch==0)
    {
      hleft.FillRandom("form1",10000);
      hright.FillRandom("form2",10000);
    }
  else
    {
      for (int i=1;i<100;++i)
	{
	  hleft.SetBinContent(i, i*(100-i)*(100-i));  // some function
	  hright.SetBinContent(i, i*i*(100-i));  // some other function
	}
    }

  TCanvas mycanvas("c1","Interp Histo",800,800);
  mycanvas.Divide(1,2);
  gStyle->SetOptStat(0);
  hleft.SetLineColor(1);
  hright.SetLineColor(2);
  hinterph.SetLineColor(4);
  hinterpv.SetLineColor(4);

  for (double x=4.5; x<5.5; x += 0.002)
    {
      mycanvas.cd(1);
      csm_interpolate_histogram(&hleft,4.5,&hright,5.5,&hinterph,x,CSM_INTERP_HORIZONTAL);
      hleft.Draw("hist");
      hright.Draw("same,hist");
      hinterph.Draw("same,hist");
      mycanvas.cd(2);
      csm_interpolate_histogram(&hleft,4.5,&hright,5.5,&hinterpv,x,CSM_INTERP_VERTICAL);
      hleft.Draw("hist");
      hright.Draw("same,hist");
      hinterpv.Draw("same,hist");
      mycanvas.Update();
      gSystem->ProcessEvents();
      gSystem->Sleep(10);
    }
}


void te2(int iswitch=0)
{
  TH2F hleft("hleft"," ",100,0,1,100,0,1);
  TH2F hright("hright"," ",100,0,1,100,0,1);
  TH2F hinterph("hinterph"," ",100,0,1,100,0,1);

  TF2 form1("form1","exp(-( (x-0.2)*(x-0.2) + (y-0.2)*(y-0.2) )/0.005)");
  TF2 form2("form2","exp(-( fabs(x-0.7) + fabs(y-0.7) )/0.05)");

  if (iswitch==0)
    {
      hleft.FillRandom("form1",100000);
      hright.FillRandom("form2",100000);
    }
  else
    {
      for (int i=1;i<=100;++i)
	{
	  for (int j=1;j<=100;++j)
	    {
	      hleft.SetBinContent(i, j, TMath::Max( 10-TMath::Abs(10-i),0)*TMath::Max( 10-TMath::Abs(10-j),0));  // some function
	      hright.SetBinContent(i, j, TMath::Max( 10-TMath::Abs(80-i),0)*TMath::Max( 10-TMath::Abs(80-j),0));  // some other function
	    }
	}
    }

  TCanvas mycanvas("c1","Interp Histo",800,800);
  gStyle->SetOptStat(0);
  hinterph.SetLineColor(4);

  for (double x=4.5; x<5.5; x += 0.002)
    {
      csm_interpolate_histogram(&hleft,4.5,&hright,5.5,&hinterph,x,CSM_INTERP_HORIZONTAL);
      hleft.Draw("cont");
      hright.Draw("same,cont");
      hinterph.Draw("same,cont");
      mycanvas.Update();
      gSystem->ProcessEvents();
      gSystem->Sleep(10);
    }
}



void tesyst(int iswitch=0)
{
  TH1F hleft("hleft"," ",100,0,1);
  TH1F hsyst("hysst"," ",100,0,1);
  TH1F hright("hright"," ",100,0,1);
  TH1F hinterph("hinterph"," ",100,0,1);
  TH1F hinterph2("hinterph2"," ",100,0,1);

  TF1 form1("form1","exp(-(x-0.2)*(x-0.2)/0.003)");
  TF1 form2("form2","exp(-(x-0.7)*(x-0.7)/0.003)");
  TF1 form3("form3","exp(-(x-0.2)*(x-0.2)/0.01)");

  if (iswitch==0)
    {
      hleft.FillRandom("form1",10000);
      hright.FillRandom("form2",10000);
      hsyst.FillRandom("form3",10000);
    }
  else
    {
      std::cout << "iswitch not implementedd yet" << std::endl;
    }

  TCanvas mycanvas("c1","Interp Histo",800,800);
  gStyle->SetOptStat(0);
  hleft.SetLineColor(1);
  hright.SetLineColor(2);
  hsyst.SetLineColor(2);
  hinterph.SetLineColor(4);
  hinterph2.SetLineColor(4);

  for (double x=4.5; x<5.1; x += 0.002)
    {
      csm_interpolate_histogram(&hleft,4.5,&hright,5.5,&hinterph,x,CSM_INTERP_HORIZONTAL);
      hleft.Draw("hist");
      hright.Draw("same,hist");
      hinterph.Draw("same,hist");
      mycanvas.Update();
      gSystem->ProcessEvents();
      gSystem->Sleep(10);
    }

  hright.SetLineStyle(2);
  for (double x=0.0; x<1.0; x += 0.002)
    {
      csm_interpolate_histogram2(&hleft,0.0,&hsyst,1.0,&hinterph,&hinterph2,x,CSM_INTERP_HORIZONTAL);
      hleft.Draw("hist");
      hsyst.Draw("same,hist");
      hright.Draw("same,hist");
      hinterph2.Draw("same,hist");
      mycanvas.Update();
      gSystem->ProcessEvents();
      gSystem->Sleep(10);
    }
  gSystem->Sleep(10000);
}
