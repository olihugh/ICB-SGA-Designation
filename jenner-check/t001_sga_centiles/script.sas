/**************************************************************************
 Bundle:  t001_sga_centiles
 Source:  "ICB SGA Designation Paper Script.sas" (data step1 + descriptive tabulations)
 Author of the SAS logic: Oliver Hugh, Perinatal Institute

 This is the centile-computation core of the paper script, run against a
 small synthetic cohort in place of the proprietary inital_db extract.
 Every growth-chart formula (FMF, Hadlock, WHO, IG21-2017, IG21-2020) and
 every SGA/LGA designation below is exactly as written in the original
 program; only the data source and the Excel output targets are swapped
 so the program runs stand-alone.
**************************************************************************/

/* Synthetic cohort standing in for db.inital_db.
   Columns match those the DATA step reads. gest is in days. */
data inital_db;
  infile datalines dsd;
  length eth $8 icb $6;
  input gest birthweight bmi parity eth $ ccb sndb icb $
        maternalheight maternalweight;
datalines;
280,3400,24.1,0,1.Brit,52.0,0.31,ICB001,164,65
259,2650,22.0,1,2.SA,8.5,-1.42,ICB001,158,55
294,3980,31.4,0,4.Bla,88.0,1.20,ICB002,170,91
238,1420,19.8,0,2.SA,2.1,-2.10,ICB002,160,51
287,3720,27.6,2,1.Brit,64.0,0.55,ICB003,167,77
266,2900,20.3,0,3.EE,26.0,-0.66,ICB003,162,53
301,4210,35.9,1,5.Oth,95.0,1.75,ICB001,171,105
245,1980,23.7,0,1.Brit,14.0,-1.05,ICB002,159,60
273,3150,25.5,1,4.Bla,41.0,-0.22,ICB003,166,70
290,3560,29.1,0,3.EE,58.0,0.20,ICB001,163,78
252,2380,21.2,0,2.SA,5.5,-1.63,ICB002,157,52
298,4050,33.0,3,1.Brit,91.0,1.34,ICB003,169,95
263,2760,18.2,0,5.Oth,18.0,-0.91,ICB001,161,47
284,3480,26.0,1,1.Brit,49.0,0.08,ICB002,165,71
277,3260,24.9,0,4.Bla,44.0,-0.14,ICB003,166,68
;
run;

data step1;
set inital_db;

if gest<24*7 then delete;
if gest>42*7 then delete;


** FMF;
gest2=gest;
bwt2=birthweight;

FMFP50 = 3.0893+0.00835*(gest2-199)-0.00002965*(gest2-199)**2 -0.00000006062*(gest2-199)**3 ;
FMFEP10 = FMFP50 - 1.282*(0.02464+0.0000564*gest2);
FMFEP3 = FMFP50 - 1.88*(0.02464+0.0000564*gest2);

FMFEP90 = FMFP50 + 1.282*(0.02464+0.0000564*gest2);
FMFEP97 = FMFP50 + 1.88*(0.02464+0.0000564*gest2);

FMFBP10 = FMFP50 - 1.282*sqrt((0.02464+0.0000564*gest2)**2 + (0.03363**2));
FMFBP90 = FMFP50 + 1.282*sqrt((0.02464+0.0000564*gest2)**2 + (0.03363**2));

if bwt2 < 10**(FMFEP10) then FMFESGA=1; else FMFESGA=0;
if bwt2 < 10**(FMFEP3) then FMFESGA3=1; else FMFESGA3=0;

if bwt2 > 10**(FMFEP90) then FMFELGA=1; else FMFELGA=0;
if bwt2 > 10**(FMFEP97) then FMFELGA97=1; else FMFELGA97=0;

  mu     = 3.0893 + 0.00835*(gest2-199) - 0.00002965*(gest2-199)**2 - 0.00000006062*(gest2-199)**3;
  sigma  = 0.02464 + 0.0000564*gest2;

 logbw  = log10(bwt2);
  z      = (logbw - mu) / sigma;
  fmfcent = round(probnorm(z) * 100, 0.1);  format fmfcent 6.2;

FMFdec = ceil(fmfcent/10);

if fmfcent=10 then fmfdec=2;
if fmfcent=90 then fmfdec=9;

** HADLOCK;

hadp50exactgenb = exp(0.578 + 0.332*(gest2/7) - 0.00354 * (gest2/7)**2);
hadp50exact280 = exp(0.578 + 0.332*(280/7) - 0.00354 * (280/7)**2);

efwth= (hadp50exactgenb /  hadp50exact280)*3619;

sndh = (bwt2- efwth )/ (0.1333* efwth);

cch = round(100 * cdf('normal', sndh), 0.1);

if cch<10 then SGAHad=1; else SGAHad=0;
if cch<3 then SGAHad3=1; else SGAHad3=0;

if cch>90 then lGAHad=1; else LGAHad=0;
if cch>97 then lGAHad97=1; else LGAHad97=0;

haddec = ceil(cch/10);

if cch=10 then haddec=2;
if cch=90 then haddec=9;

** WHO;


gest2_2=(gest2)**2;
gest2_3=(gest2)**3;


WHOfetp_25bwt=exp(9.65473E-08*gest2_3-0.000128714*gest2_2+0.057725657*gest2-0.262514302);
WHOfetp5bwt=exp(7.05336E-08*gest2_3-0.000115928*gest2_2+0.055813478*gest2-0.133011937);
WHOfetp10bwt=exp(7.41328E-08*gest2_3-0.000115165*gest2_2+0.055185419*gest2-0.016656143);
WHOfetp25bwt=exp(8.38185E-08*gest2_3-0.000120657*gest2_2+0.056213613*gest2-0.008739158);
WHOfetp50bwt=exp(7.43779E-08*gest2_3-0.00011459*gest2_2+0.055024219*gest2+0.138046844);
WHOfetp75bwt=exp(5.99688E-08*gest2_3-0.000105783*gest2_2+0.053299317*gest2+0.319301605);
WHOfetp90bwt=exp(6.47403E-08*gest2_3-0.000109177*gest2_2+0.05402464*gest2+0.33895772);
WHOfetp95bwt=exp(7.89464E-08*gest2_3-0.000117444*gest2_2+0.055492027*gest2+0.305653135);
WHOfetp_975bwt=exp(7.34664E-08*gest2_3-0.000113307*gest2_2+0.054594203*gest2+0.399670065);



if bwt2<WHOfetp5bwt then whocent=5-2.5*(WHOfetp5bwt-bwt2)/(WHOfetp5bwt-WHOfetp_25bwt);
if bwt2>=WHOfetp5bwt and bwt2<WHOfetp10bwt then whocent=10-5*(WHOfetp10bwt-bwt2)/(WHOfetp10bwt-WHOfetp5bwt);
if bwt2>=WHOfetp10bwt and bwt2<WHOfetp25bwt then whocent=25-15*(WHOfetp25bwt-bwt2)/(WHOfetp25bwt-WHOfetp10bwt);
if bwt2>=WHOfetp25bwt and bwt2<WHOfetp50bwt then whocent=50-25*(WHOfetp50bwt-bwt2)/(WHOfetp50bwt-WHOfetp25bwt);
if bwt2>=WHOfetp50bwt and bwt2<WHOfetp75bwt then whocent=75-25*(WHOfetp75bwt-bwt2)/(WHOfetp75bwt-WHOfetp50bwt);
if bwt2>=WHOfetp75bwt and bwt2<WHOfetp90bwt then whocent=90-15*(WHOfetp90bwt-bwt2)/(WHOfetp90bwt-WHOfetp75bwt);
if bwt2>=WHOfetp90bwt and bwt2<WHOfetp95bwt then whocent=95-5*(WHOfetp95bwt-bwt2)/(WHOfetp95bwt-WHOfetp90bwt);
if bwt2>=WHOfetp95bwt  then whocent=97.5-2.5*(WHOfetp_975bwt-bwt2)/(WHOfetp_975bwt-WHOfetp95bwt);
if whocent<0 then whocent=0; if whocent>100 then whocent=100;
whocent=round(whocent,0.1);

if bwt2<WHOfetp10bwt then WHOSGA=1; else WHOSGA=0;
if whocent<3 then WHOSGA3=1; else WHOSGA3=0;
if whocent>90 then WHOlGA=1; else WHOlGA=0;
if whocent>97 then WHOlGA97=1; else WHOlGA97=0;

whodec = ceil(whocent/10);

if whocent=10 then whodec=2;
if whocent=90 then whodec=9;

** GROW Lite;

hadp50exactgenb = exp(0.578 + 0.332*(gest2/7) - 0.00354 * (gest2/7)**2);
hadp50exactmed = exp(0.578 + 0.332*(280/7) - 0.00354 * (280/7)**2);

ccnhs = round(100 * cdf('normal', sndb), 0.1);
if ccnhs<10 then nhssga=1; else nhssga=0;
if ccnhs<3 then nhssga3=1; else nhssga3=0;

if ccnhs>90 then nhslga=1; else nhslga=0;
if ccnhs>97 then nhslga97=1; else nhslga97=0;

nhsdec = ceil(ccnhs/10);

if ccnhs=10 then dec=2;
if ccnhs=90 then dec=9;

** GROW, already recorded as 'ccb' at time of care;

if ccb<10 then sga=1; else sga=0;
if ccb<3 then sga3=1; else sga3=0;

if ccb>90 then lga=1; else lga=0;
if ccb>97 then lga97=1; else lga97=0;


dec = ceil(ccb/10);

if ccb=10 then dec=2;
if ccb=90 then dec=9;


** IG21 2017, modelled from tables included within paper;

gest2_2=gest2**2;
gest2_3=gest2**3;
gest2_4=gest2**4;

IGfetp10bwt=exp(4.628946519+0.003180572*gest2+3.96909E-05*gest2_2+1.09034E-07*gest2_3+-5.02005E-10*gest2_4);
IGfetp90bwt=exp(6.30147567+-0.029993015*gest2+0.000281046*gest2_2+-5.59424E-07*gest2_3+9.82821E-11*gest2_4);
IGfetp50bwt=exp(5.449988096+-0.011975144*gest2+0.000141117*gest2_2+-1.38218E-07*gest2_3+-3.27535E-10*gest2_4);
IGfetp3bwt=exp(4.79790337-0.000977031*gest2+7.67137E-05*gest2_2-5.02226E-08*gest2_3-2.58745E-10*gest2_4)			 ;
IGfetp97bwt=exp(6.786773085+-0.040814907*gest2+0.000369575*gest2_2+-8.47234E-07*gest2_3+4.21186E-10*gest2_4);

if bwt2<IGfetp10bwt then SGAIGFet=1 ; else SGAIGFet=0;
if bwt2<IGfetp3bwt then SGAIGFet3=1 ; else SGAIGFet3=0;

if bwt2>IGfetp90bwt then LGAIGFet=1 ; else LGAIGFet=0;
if bwt2>IGfetp97bwt then LGAIGFet97=1 ; else LGAIGFet97=0;

  lnP03 = log(IGfetp3bwt);
  lnP10 = log(IGfetp10bwt);
  lnP50 = log(IGfetp50bwt);
  lnP90 = log(IGfetp90bwt);
  lnP97 = log(IGfetp97bwt);

   P03 = exp(4.79790337 + -0.000977031*gest2 + 7.67137E-05*gest2_2 - 5.02226E-08*gest2_3 - 2.58745E-10*gest2_4);
  P10 = exp(4.628946519 + 0.003180572*gest2 + 3.96909E-05*gest2_2 + 1.09034E-07*gest2_3 + -5.02005E-10*gest2_4);
  P50 = exp(5.449988096 + -0.011975144*gest2 + 0.000141117*gest2_2 + -1.38218E-07*gest2_3 + -3.27535E-10*gest2_4);
  P90 = exp(6.30147567  + -0.029993015*gest2 + 0.000281046*gest2_2 + -5.59424E-07*gest2_3 + 9.82821E-11*gest2_4);
  P97 = exp(6.786773085 + -0.040814907*gest2 + 0.000369575*gest2_2 + -8.47234E-07*gest2_3 + 4.21186E-10*gest2_4);

  /* Log scale */
  lnP03 = log(P03); lnP10 = log(P10); lnP50 = log(P50); lnP90 = log(P90); lnP97 = log(P97);
  lnbw  = log(bwt2) ;

  length band $8;
  igcent = .;

  if missing(lnbw) then do;
    igcent = .;
  end;
  else do;
    /* Piecewise linear interpolation on ln-scale */
    if lnbw <= lnP03 then do;               /* extrapolate below 3rd using 3-10 slope */
      igcent = 3 + (lnbw - lnP03) * (7 / (lnP10 - lnP03));
      band = '<3';
    end;
    else if lnbw <= lnP10 then do;          /* 3 - 10 (7 pp) */
      igcent = 3 + 7  * (lnbw - lnP03) / (lnP10 - lnP03);
      band = '3-10';
    end;
    else if lnbw <= lnP50 then do;          /* 10 - 50 (40 pp) */
      igcent = 10 + 40 * (lnbw - lnP10) / (lnP50 - lnP10);
      band = '10-50';
    end;
    else if lnbw <= lnP90 then do;          /* 50 - 90 (40 pp) */
      igcent = 50 + 40 * (lnbw - lnP50) / (lnP90 - lnP50);
      band = '50-90';
    end;
    else if lnbw <= lnP97 then do;          /* 90 - 97 (7 pp) */
      igcent = 90 + 7  * (lnbw - lnP90) / (lnP97 - lnP90);
      band = '90-97';
    end;
    else do;                                /* extrapolate above 97th using 90-97 slope */
      igcent = 97 + (lnbw - lnP97) * (3 / (lnP97 - lnP90));
      band = '>97';
    end;

    /* tidy bounds and exact ties */
    if abs(lnbw - lnP10) < 1e-12 then igcent = 10;
    if abs(lnbw - lnP50) < 1e-12 then igcent = 50;
    igcent = max(0, min(100, igcent));
  end;

  format igcent 6.2;

igcent = round(igcent,0.1);
igdec = ceil(igcent/10);

if igcent=10 then igdec=2;
if igcent=90 then igdec=9;

** IG21 2020;

gestw=gest2/7;
mean = -2.42272+1.86478*gestw**0.5 - 0.0000193299*gestw**3;
skew = 9.43643+9.41579*(gestw/10)**-2 - 83.5422*log(gestw/10)*(gestw/10)**-2;
cv = 0.0193557 + 0.0310716*(gestw/10)**-2 - 0.0657587*log(gestw/10)*(gestw/10)**-2;

z_score = ((cv*skew)**-1) * ((log(bwt2)/mean)**skew -1 );
ig_cent= round(probnorm(z_score)*100,0.1);

if ig_cent<10 then sgaigfet20=1; else sgaigfet20=0;
if ig_cent<3 then sgaigfet203=1; else sgaigfet203=0;

if ig_cent>90 then lgaigfet20=1; else lgaigfet20=0;
if ig_cent>97 then lgaigfet2097=1; else lgaigfet2097=0;

ig_dec = ceil(ig_cent/10);

if ig_cent=10 then ig_dec=2;
if ig_cent=90 then ig_dec=9;

*** Coding of descriptives;

if bmi<18.5 then bmi185=1; else bmi185=0;
if bmi>30 then bmi30=1; else bmi30=0;
if bmi>35 then bmi35=1; else bmi35=0;

if gest<37*7 then prem37=1; else prem37=0;
if gest<34*7 then prem34=1; else prem34=0;


if parity=0 then nullip=1; else nullip=0;

if eth="1.Brit" then brit=1; else brit=0;
if eth="2.SA" then SA=1; else SA=0;
if eth="3.EE" then EE=1; else EE=0;
if eth="4.Bla" then bla=1; else bla=0;
if eth="5.Oth" then oth=1; else oth=0;

if sgahad=1 and sga=0 then hadonly=1; else hadonly=0;
if sgahad=1 and sga=1 then both=1; else both=0;
if sgahad=0 and sga=1 then custonly=1; else custonly=0;
if sgahad=0 and sga=0 then neither=1; else neither=0;

if haddec=0 then haddec=1;
if igdec=0 then igdec=1;
if ig_dec=0 then ig_dec=1;
if whodec=0 then whodec=1;
if fmfdec=0 then fmfdec=1;
if dec=0 then dec=1;
if nhsdec=0 then nhsdec=1;

run;

** Per-birth centiles produced by each growth chart;
proc print data=step1 noobs;
var icb gest birthweight fmfcent cch whocent igcent ig_cent ccb;
run;

** SGA/LGA designation counts across the seven charts;
proc tabulate data=step1;
class icb;
var SGAHad SGAIGFet sgaigfet20 WHOSGA FMFESGA sga nhssga;
table icb='', n sum*(SGAHad SGAIGFet sgaigfet20 WHOSGA FMFESGA sga nhssga)*f=6.0;
run;

** Decile agreement summary between charts;
proc means data=step1 n mean min max maxdec=1;
var fmfcent cch whocent igcent ig_cent;
run;
