/**************************************************************************
 Program:    Paper 39.sas
 Author:     Oliver Hugh
 Organisation: Perinatal Institute
 Contact:    ohugh@perinatal.org.uk
 Created:    07-Oct-2025
 Updated:    07-Oct-2025
 SAS Version: 9.4
 Platform:   Windows 

 Description:
   This program is used to produce the tables and figures referenced in  "Designation of small-for-gestational age according to 7 fetal growth charts in 
England’s National Health Service: retrospective analysis of 3.2 million births" by Gardosi et al, 2026. 


 Inputs:
   - inital_db.sas7bdat: dataset containing 3.2 million births with identification of ICB (available on reasonable request)
 Outputs:
	- Figure 1.xlsx
	- OTHER FIGURES AND TABLES.xlsx

 License:
   © 2026 Perinatal Institute. All rights reserved.
   This program is provided for research and educational use only.
   Redistribution or modification requires prior written permission.

 Change Log:
   Date        Author        Description
   ----------  ------------- ------------------------------------------
   07Oct2025   O. Hugh       Initial version
**************************************************************************/


libname db "P:\Large Databases\SAS\Paper 39\01oct25";

data step1;
set db.inital_db;

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
    if lnbw <= lnP03 then do;               /* extrapolate below 3rd using 3–10 slope */
      igcent = 3 + (lnbw - lnP03) * (7 / (lnP10 - lnP03));
      band = '<3';
    end;
    else if lnbw <= lnP10 then do;          /* 3 ? 10 (7 pp) */
      igcent = 3 + 7  * (lnbw - lnP03) / (lnP10 - lnP03);
      band = '3-10';
    end;
    else if lnbw <= lnP50 then do;          /* 10 ? 50 (40 pp) */
      igcent = 10 + 40 * (lnbw - lnP10) / (lnP50 - lnP10);
      band = '10-50';
    end;
    else if lnbw <= lnP90 then do;          /* 50 ? 90 (40 pp) */
      igcent = 50 + 40 * (lnbw - lnP50) / (lnP90 - lnP50);
      band = '50-90';
    end;
    else if lnbw <= lnP97 then do;          /* 90 ? 97 (7 pp) */
      igcent = 90 + 7  * (lnbw - lnP90) / (lnP97 - lnP90);
      band = '90-97';
    end;
    else do;                                /* extrapolate above 97th using 90–97 slope */
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

** GROUPING BY ICB;

proc tabulate data=step1 out=grouped;
class haddec igdec ig_dec whodec fmfdec dec nhsdec  icb;
table  icb='', n (haddec igdec ig_dec whodec fmfdec dec nhsdec)*n;
run;

data overall;
set grouped;
where _TYPE_="00000001";
run;

proc sql noprint;
create table Final_DB as
select a.*, b.n as tot
from overall a, grouped b 
where a.icb=b.icb and a._TYPE_^="00000001" and b.n>=1000;
quit;

data Final_DB;
set Final_DB;
perc = n/tot*100;
run;


ods excel file="P:\Oliver folder\SAS\Output\&sysdate. Figure 1.xlsx" style=styles.htmlblue
options(sheet_name="Hadlock" frozen_headers="1") ;

proc tabulate data=Final_DB;
where haddec^=.;
var perc;
class haddec;
table haddec='', perc*mean*f=6.1;
run;

ods excel options(sheet_name="IG 17" frozen_headers="1") ;

proc tabulate data=Final_DB;
where igdec^=.;
var perc;
class igdec;
table igdec='', perc*mean*f=6.1;
run;

ods excel options(sheet_name="IG 20" frozen_headers="1") ;


proc tabulate data=Final_DB;
where ig_dec^=.;
var perc;
class ig_dec;
table ig_dec='', perc*mean*f=6.1;
run;

ods excel options(sheet_name="WHO" frozen_headers="1") ;


proc tabulate data=Final_DB;
where whodec^=.;
var perc;
class whodec;
table whodec='', perc*mean*f=6.1;
run;

ods excel options(sheet_name="FMF" frozen_headers="1") ;


proc tabulate data=Final_DB;
where fmfdec^=.;
var perc;
class fmfdec;
table fmfdec='', perc*mean*f=6.1;
run;

ods excel options(sheet_name="GROW" frozen_headers="1") ;


proc tabulate data=Final_DB;
where dec^=.;
var perc;
class dec;
table dec='', perc*mean*f=6.1;
run;

ods excel options(sheet_name="GROWL" frozen_headers="1") ;


proc tabulate data=Final_DB;
where nhsdec^=.;
var perc;
class nhsdec;
table nhsdec='', perc*mean*f=6.1;
run;

ods excel close;

** GROUPING BY ICB**;
proc sql noprint;
create table step2 as
select distinct icb, n(icb) as n, mean(bmi) as bmi_mean, median(bmi) as bmi_p50, mean(maternalheight) as height_mean, mean(maternalweight) as weight_mean, median(maternalweight) as weight_p50, 
sum(bmi185) as bmi185, sum(bmi30) as bmi30, sum(bmi35) as bmi35, sum(prem37) as prem37, sum(prem34) as prem34,
sum(brit) as brit, sum(SA) as sa, sum(ee) as ee, sum(bla) as bla, sum(oth) as oth, sum(nullip) as nullip
from step1 
group by icb;
quit;

data step2;
set step2;

bmi185_ = bmi185/n*100;
bmi30_ = bmi30/n*100;
bmi35_ = bmi35/n*100;

prem37_ = prem37/n*100;
prem34_ = prem34/n*100;

Brit_ = brit/n*100;
sa_ = sa/n*100;
ee_ = ee/n*100;
bla_ = bla/n*100;
oth_ = oth/n*100;

nullip_ = nullip/n*100;

run;

ods excel file="P:\Oliver folder\SAS\Output\&sysdate. OTHER FIGURES AND TABLES.xlsx" style=styles.htmlblue
options(sheet_name="Descriptives" frozen_headers="1") ;

proc tabulate data=step2;
where n>=1000;
var Brit_ sa_ ee_ bla_ oth_ nullip_ bmi_mean bmi_p50 bmi185_ bmi30_ bmi35_ height_mean weight_mean weight_p50 prem37_ ;
table Brit_ sa_ ee_ bla_ oth_ nullip_ height_mean weight_mean weight_p50 bmi_mean bmi_p50 bmi185_ bmi30_ bmi35_  prem37_ , min max;
run;

proc tabulate data=step1;
class  parity bmi185 bmi30 bmi35 prem37 prem34 brit ee sa bla oth;
var maternalheight maternalweight bmi gest;
table  parity  bmi185 bmi30 bmi35 prem37 prem34 brit ee sa bla oth, n pctn;
table maternalheight maternalweight bmi gest, mean std p50 p25 p75;
run;

proc tabulate data=STEP1;
class nhssga WHOSGA  SGAHad FMFESGA FMFBSGA sga ;
table FMFESGA FMFBSGA WHOSGA     SGAHad  sga nhssga all, n;
run;

proc tabulate data=STEP1;
class nhssga WHOSGA  SGAHad FMFESGA FMFBSGA sga eth;
table FMFESGA FMFBSGA WHOSGA     SGAHad  sga nhssga all, n*eth;
run;

proc print data=STEP1 noobs;
var icb n brit_ sa_ ee_ bla_ oth_ nullip_ bmi_mean bmi_p50 height_mean weight_mean weight_p50 bmi185_ bmi30_ bmi35_ prem37_;
run;

ods excel options(sheet_name="SGA" frozen_headers="1") ;


proc tabulate data=STEP1;
class icb;
var SGAHad SGAIGFet sgaigfet20  WHOSGA FMFESGA sga nhssga  brit bmi bmi30 bmi35 ;
table icb='', n (p50 mean)*bmi sum*( bmi30 bmi35 brit SGAHad SGAIGFet sgaigfet20 WHOSGA FMFESGA sga nhssga )*f=8.0;
run;

ods excel options(sheet_name="SGA3" frozen_headers="1") ;


proc tabulate data=STEP1;
class icb;
var SGAHad3 SGAIGFet3 sgaigfet203 WHOSGA3 FMFESGA3 sga3 nhssga3  brit bmi bmi30 bmi35;
table icb='', n (p50 mean)*bmi sum*(bmi30 bmi35 brit SGAHad3 SGAIGFet3 sgaigfet203 WHOSGA3 FMFESGA3 sga3 nhssga3 )*f=8.0;
run;

ods excel options(sheet_name="LGA" frozen_headers="1") ;


proc tabulate data=STEP1;
class icb;
var LGAHad LGAIGFet Lgaigfet20  WHOLGA FMFELGA Lga nhsLga  brit bmi bmi30 bmi35;
table icb='', n (p50 mean)*bmi sum*(bmi30 bmi35 brit LGAHad LGAIGFet Lgaigfet20 WHOLGA FMFELGA Lga nhsLga )*f=8.0;
run;

ods excel options(sheet_name="LGA97" frozen_headers="1") ;


proc tabulate data=STEP1;
class icb;
var LGAHad97 LGAIGFet97 Lgaigfet2097  WHOLGA97 FMFELGA97 Lga97 nhsLga97  brit bmi bmi30 bmi35;
table icb='', n (p50 mean)*bmi sum*(bmi30 bmi35 brit LGAHad97 LGAIGFet97 Lgaigfet2097  WHOLGA97 FMFELGA97 Lga97 nhsLga97 )*f=8.0;
run;

ods excel options(sheet_name="SGATERM" frozen_headers="1") ;


proc tabulate data=Step1;
where gest>=37*7;
class icb;
var SGAHad SGAIGFet sgaigfet20  WHOSGA FMFESGA sga nhssga  brit bmi bmi30 bmi35;
table icb='', n (p50 mean)*bmi sum*(bmi30 bmi35 brit SGAHad SGAIGFet sgaigfet20 WHOSGA FMFESGA sga nhssga )*f=8.0;
run;

ods excel options(sheet_name="SGA3TERM" frozen_headers="1") ;


proc tabulate data=Step1;
where gest>=37*7;
class icb;
var SGAHad3 SGAIGFet3 sgaigfet203 WHOSGA3 FMFESGA3 sga3 nhssga3  brit bmi bmi30 bmi35;
table icb='', n (p50 mean)*bmi sum*(bmi30 bmi35 brit SGAHad3 SGAIGFet3 sgaigfet203 WHOSGA3 FMFESGA3 sga3 nhssga3 )*f=8.0;
run;


ods excel close;

**** Linear Regression - Table 4; 

data raw; * adjusting ICBs to be on similar scale; 
set step2;
eng = __English/10;
hei = Mean_height*2;
wei= Mean_weight;
run;

ods graphics off;

proc reg data=raw;
   model Hadlock = wei / clb  ;
run;
quit;

proc reg data=raw;
   model IG21_2017 = wei / clb  ;
run;
quit;

proc reg data=raw;
   model IG21_2020 = wei / clb  ;
run;
quit;

proc reg data=raw;
   model who = wei / clb  ;
run;
quit;

proc reg data=raw;
   model fmf = wei / clb  ;
run;
quit;

proc reg data=raw;
   model grow = wei / clb  ;
run;
quit;

proc reg data=raw;
   model grow_lite = wei / clb  ;
run;
quit;
