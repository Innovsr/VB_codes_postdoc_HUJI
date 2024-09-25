!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                            Chemical Insight Structure Sets                             !!
!!                                                                                        !!
!!                               Written by: Dr. Sourav Roy                               !!
!!                                                                                        !!
!!                    School Of Pharmacy, Hebrew University of Jerusalem                  !!
!!                                                                                        !!
!!                                    Jerusalem, Israel                                   !!
!!                                                                                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module commondat1
implicit none

integer::orbs(50,500),orbsl,freezorb(100),frzn,orbn(50),sl1(10),orbsym(20,20),norbsym(50),frag_at(20,20)&
,frag_atn(20),atm_nb_orbs1(20),atm_nb_orbs(20),ovlp_int,nbasatom(20),nbasis,ovlp_bas_num(20,10),lnum&
,num_act_orb_typ(20),ind_mat(1000,1000)
character(len=5)::symat(20),all_atbas(500)
real*8::symatno(20)
double precision::orb_ovlp_mat(20,20),orb_ovlp_mat1(20,20),ovlp_mat_norm(5000,100)

save
end module commondat1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module commondat
implicit none

integer::nao,nae,niao,mult,nlast,nlp,tncs,tnis,is,atom,nabsym,ndet,symm,repl,itb,syb,nnb,rad &
,flgst,noqult,rplc,flg_ion,flg_cov,imbd,nfset,noq0,noq1,noq2,noq3,noq,niabd,nialp,niach,serial &
,nrad,choice,pbond(1000),mbond(1000),nnnatom,nnat_bond(100,2),prio_rad(20,100),norad(20),prad &
,nlpset,lp(100),plpair(100,10),flg1,uoptstr,mnbond,nmbond,main_bond(100),iabd(50),ialp(50),iach(50)
integer::atoset(200,20),at_num(88),at_covrad(88),key_frag,frag_cntr(10)&
,at_ab_symset(20,20),tot_ndet,nrs,mset,loopsymsc(15000) &
,atn(50),absyn(50),qult(100),radical,actv_atom(30),num_frag_cntr,input_flg,ovopt,vpt,all_at_num(100)&
,strt_struc(100,15),nstrt,repl_struc(100,15),nnmat_act(100,2),nnmat_inact(100,2),nactorb,active_orbs(20),&
active_atoms(30),nnat_bond_inact(100,2),nnatominact,tot_atom,biasmat(100),sig_sym_flg,valence_state(54)&
,bond_count(1000,100),totrum,vacorb
character(len=2),public ::at_list(88),at_list_bold(88),int_num(40)
integer,public::atm_nb_sym(20),nalpha,nbeta,str_det_sec(15000,1000),num_norbsym1(100),qflg
!real,dimension(:,:),allocatable::ovlp_mat
real*8::ovval,dist_rel_mat(20,20),biasval(20)&
,dist_nnat(20)
real*8::prime_num(142),dist_mat(20,20),iab_length
character(len=6),public::symtype,norbsym1(100,20)
double precision::dist_act_rel_mat(20,20)
integer::strdet(10000),detmnt(10000,15),det_sign(10000),Rumwrite,set_order

common /info/nlp,tncs,tnis,is,atom,nabsym,atn,absyn,active_atoms,nao,nae,niao,mult,nlast,set_order,loopsymsc
common /mat/dist_mat,flgst,noqult,rplc,nactorb,nstrt,serial,active_orbs,atm_nb_sym,dist_rel_mat,str_det_sec,biasval

common /stct/mnbond,nmbond,nfset,at_num,itb,syb,nnb,rad,prad,noq,noq1,noq2,noq3,dist_act_rel_mat&
,atoset,at_ab_symset,repl,radical,tot_atom,input_flg,biasmat,dist_nnat,Rumwrite,qflg

common /srr/key_frag,num_frag_cntr,imbd,qult,nnat_bond,nnnatom,frag_cntr
common /str1/all_at_num,sig_sym_flg,uoptstr,repl_struc,ovval,actv_atom,noq0,nlpset,lp,plpair,strt_struc,flg1,ndet&
,ovopt,vpt,symtype

DATA valence_state/1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,&
4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5/

DATA at_list/'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al',&
'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni',&
'Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc',&
'Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce',&
'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',&
'W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra'/

DATA at_list_bold/'H','HE','LI','BE','B','C','N','O','F','NE','NA','MG','AL',&
'SL','P','S','CL','AR','K','CA','SC','TI','V','CR','MN','FE','CO','NR',&
'CU','ZN','GA','GE','AS','SE','BR','KR','RB','SR','Y','ZR','NB','MO','TC',&
'RU','RH','PD','AG','CD','IN','SN','SB','TE','I','XE','CS','BA','LA','CE',&
'PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA',&
'W','RE','OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN','FR','RA'/

DATA int_num/'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16',&
'17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33',&
'34','35','36','37','38','39','40'/

DATA prime_num/11.0, 13.0, 17.0, 19.0, 23.0, 29.0, 31.0, 37.0, 41.0, 43.0, 47.0, 53.0, 59.0, 61.0&
, 67.0, 71.0, 73.0, 79.0, 83.0, 89.0, 97.0, 101.0, 103.0, 107.0, 109.0, 113.0, 127.0, 131.0, 137.0, 139.0,&
 149.0, 151.0, 157.0, 163.0, 167.0, 173.0, 179.0, 181.0, 191.0, 193.0, 197.0, 199.0, 211.0, 223.0, 227.0, 229.0,&
 233.0, 239.0, 241.0, 251.0, 257.0, 263.0, 269.0, 271.0, 277.0, 281.0, 283.0, 293.0, 307.0, 311.0, 313.0, 317.0,&
 331.0, 337.0, 347.0, 349.0, 353.0, 359.0, 367.0, 373.0, 379.0, 383.0, 389.0, 397.0, 401.0, 409.0, 419.0, 421.0,&
 432.0, 433.0, 439.0, 443.0, 449.0, 457.0, 461.0, 463.0, 467.0, 479.0, 487.0, 491.0, 499.0, 503.0, 509.0, 521.0,&
 523.0, 541.0, 547.0, 557.0, 563.0, 569.0, 571.0, 577.0, 587.0, 593.0, 599.0, 601.0, 607.0, 613.0, 617.0, 619.0,&
 631.0, 641.0, 643.0, 647.0, 653.0, 659.0, 661.0, 673.0, 677.0, 683.0, 691.0, 701.0, 709.0, 719.0, 727.0, 733.0,&
 739.0, 743.0, 751.0, 757.0, 761.0, 769.0, 773.0, 787.0, 797.0, 809.0, 811.0, 821.0, 823.0, 827.0, 829.0, 839.0/

!!!! Covalent radii (in pm unit)  taken from 'Webelement' ... from the link
!"http://crystalmaker.com/support/tutorials/atomic-radii/index.html"

DATA at_covrad /37, 32, 134, 90, 82, 77, 75, 73, 71, 69, 154, 130, 118, 111,&
106, 102, 99, 97, 196, 174, 144, 136, 125, 127, 139, 125, 126, 121, 138, 131,&
126, 122, 119, 116, 114, 110, 211, 192, 162, 148, 137, 145, 156, 126, 135, 131,153,&
148, 144, 141, 138, 135, 133, 130, 225, 198, 169, 204, 203, 201, 199, 198, 198,&
196, 194, 192, 192, 189, 190, 187, 160, 150, 138, 146, 159, 128, 137,&
128, 144, 149, 148, 147, 146, 140, 150, 145, 260, 221/
save

end module commondat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use commondat
use commondat1
implicit none

integer::wig1,nnae,STDOUT,argnum,nstr7,str5(2000,20)
logical :: fileexists
character(len=35)::inputfilename
common/str/str5,nstr7

open(unit=7,file='structures.dat',status='unknown')
!open(unit=9,file='inputfile.xmi',status='unknown')
!open(unit=9,file='structure_set.dat',status='unknown')
open(unit=10,file='structure_set_1.dat',status='unknown')
open(unit=23,file='Rumer_Sets_all.dat',status='unknown')
open(unit=35,file='quality_str.dat',status='unknown')

serial=0
input_flg=0
nstr7=0
Rumwrite=0

argnum=iargc()
if(argnum.eq.0)then
inputfilename='input.dat'
input_flg=1
INQUIRE(FILE=TRIM(inputfilename),EXIST=fileexists)
IF (fileexists) THEN
print*,'Reading from input.dat file.'
call input
ELSE
PRINT*,'SORRY input.dat file does not exist and you also did not put .xmi file'
stop
ENDIF
endif

if(input_flg.eq.0)call read_xmi


if(flgst.eq.1.or.flgst.eq.2) then
call cov_struc
endif
if(flgst.eq.1.or.flgst.eq.3) then
call ion_struc
endif

call close_file
stop
end program main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_xmi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

integer::argnum,l,MDP,sl(50),k,ll,j,i,i1,io,endflg(50),l1,l2,gap(10,100),keyn(15),least&
,l3,l4,k1,bfimat(5,20),gap1(20),dash(20),bfiln(5),&
lastlbfi,ufrzorbn,orbn1,orbn2,gap2(200),k2,k3,k4,k5,k6,k7,nfrag
integer::chinst,orbsl1(50),k8,k9,sl2(10),MDP1
logical :: fileexists
character(len=35)::inputfilename,keywd(15,50),line7(50),line
character(len=4)::lname
character(len=55)::sttr,charst(500),string(4,50)
character(len=90)::line3,line2,line4,line5,line6
character(len=200)::lowercase
integer::atsymset(20,20),nsym,syn(50),at_sym(50),fragorb,ifrag,iorb,ibfi

common/ats/atsymset,nsym,syn,at_sym
print*,'enter read_xmi'

noq0=100
noq1=100
noq2=100
noq3=100
symm=0
ovval=1.000
ovopt=0
nlpset=0
nfset=0
flg1=1
mult=1




flgst=1
key_frag=0


!if(argnum.eq.1)then
!endif

if(input_flg.eq.1)goto 500

call getarg(1,inputfilename)
l=len(TRIM(inputfilename))
lname=inputfilename(l-3:l)

if(lname.ne.'.xmi') then
print*,'please put the .xmi file as input'
stop
endif
INQUIRE(FILE=TRIM(inputfilename),EXIST=fileexists)
IF (fileexists) THEN
open(unit=21,file=TRIM(inputfilename),status='old')
ELSE
PRINT*,'SORRY This input file does not exist or you may not provide the &
filename at all'
stop
ENDIF
write(9,*)inputfilename
write(23,*)inputfilename

MDP=0
do
read(21,'(a)',iostat=io)

if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(21)
!print*,'MDP',MDP

do i=1,MDP-1
read(21,'(a)')charst(i)
if(index(charst(i),'$end').ne.0)then
MDP1=i
endif
!!!write(9,*)charst(i)
!write(23,*)charst(i)
!sttr(i)=charst(i)(1:4)
!print*,charst(i)
enddo

rewind(21)

do i=1,10
sl1(i)=0
enddo
l=0
l1=0

do i=1,MDP1

if(index(charst(i),'#').eq.1.or.index(charst(i),';').eq.1)goto 789
if(index(charst(i),'#').gt.1)then
ll=index(charst(i),'#')
do k=1,ll
line=charst(i)(k:k)
if(line.eq.'')then
!print*,'charst(i)(1:k-1)',charst(i)(1:k-1)
if(trim(charst(i)(1:k-1)).ne.'')then
charst(i)=charst(i)(1:ll-1)
goto 111
endif
endif
enddo
endif

if(index(charst(i),';').gt.1)then
ll=index(charst(i),';')
do k=1,ll
line=charst(i)(k:k)
if(line.eq.'')then
!print*,'charst(i)(1:k-1)',charst(i)(1:k-1)
if(trim(charst(i)(1:k-1)).ne.'')then
charst(i)=charst(i)(1:ll-1)
goto 111
endif
endif
enddo
endif


111 charst(i)=lowercase(charst(i))

!ll=len(trim(charst(i)))
!do k=1,ll+1
!line5=charst(i)(k:k)
!if(line5.eq.'')then
!sttr=trim(charst(i)(1:k-1))
!if(sttr.eq.'$ctrl'.or.sttr.eq.'$end'.or.sttr.eq.'$orb'.or.sttr.eq.'$frag'.or.sttr.eq.'$bfi'&
!.or.sttr.eq.'$CTRL'.or.sttr.eq.'$END'.or.sttr.eq.'$ORB'.or.sttr.eq.'$FRAG'.or.sttr.eq.'$BFI') then
if(index(charst(i),'$ctrl').ne.0.or.index(charst(i),'$end').ne.0&
.or.index(charst(i),'$orb').ne.0.or.index(charst(i),'$frag').ne.0&
.or.index(charst(i),'$bfi').ne.0) then
l=l+1
if(index(charst(i),'$end').ne.0)then
!if(sttr.eq.'$end'.or.sttr.eq.'$END')then
endflg(l)=1
else
endflg(l)=0
endif
if(endflg(l).ne.endflg(l-1).or.l.eq.1)then
l1=l1+1
sl(l1)=i
!if(sttr.eq.'$ctrl')sl1(l1)=1
if(index(charst(i),'$ctrl').ne.0)sl1(l1)=1
!if(sttr.eq.'$CTRL')sl1(l1)=1
!if(sttr.eq.'$frag')sl1(l1)=2
if(index(charst(i),'$frag').ne.0)sl1(l1)=2
!if(sttr.eq.'$FRAG')sl1(l1)=2
!if(sttr.eq.'$orb')sl1(l1)=3
if(index(charst(i),'$orb').ne.0)sl1(l1)=3
!if(sttr.eq.'$ORB')sl1(l1)=3
!if(sttr.eq.'$bfi')sl1(l1)=4
if(index(charst(i),'$bfi').ne.0)sl1(l1)=4
!if(sttr.eq.'$BFI')sl1(l1)=4
endif
if(endflg(l).eq.0.and.endflg(l-1).eq.0)then
print*,'You might forget to put $end somewhere in the .xmi file'
stop
endif
endif
!endif
!enddo
789 enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! if fragment section present then orbital section should be read before
!!!!!!!  fragment section

fragorb=0
ifrag=2
iorb=3
k1=0
k2=0
do i=1,10
sl2(i)=0
enddo
do i=1,l1
if(sl1(i).eq.2)then
k1=1
k3=i
endif
enddo
do i=1,l1
if(sl1(i).eq.3)then
k2=1
k4=i
endif
enddo
fragorb=k1*k2
if(fragorb.eq.1)then
sl2(k3)=sl(k4)
sl2(k3+1)=sl(k4+1)
sl2(k4)=sl(k3)
sl2(k4+1)=sl(k3+1)
ifrag=3
iorb=2
do i=1,l1
if(sl2(i).ne.0)then
sl(i)=sl2(i)
endif
enddo
endif

!do i=1,l1
!print*,'sl(i)',i,sl(i)
!enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
l=0
l4=0
do i=1,l1,2
l4=l4+1
l3=0
do j=sl(i)+1,sl(i+1)-1
if(sl1(i).eq.1)then
l=l+1
l2=1
gap(1,1)=0
ll=len(trim(charst(j)))
do k=1,ll+1
line5=charst(j)(k:k)
if(line5.eq.'#')goto 110
if(line5.eq.'')then
l2=l2+1
gap(l,l2)=k
endif
enddo
110 do k=1,l2
if((gap(l,k+1)-gap(l,k)).ne.1)then
l3=l3+1
keywd(l4,l3)=trim(charst(j)(gap(l,k):gap(l,k+1)-1))
endif
enddo
endif
enddo

keyn(l4)=l3
do k=1,l3
!if(trim(keywd(l4,k)(2:7)).eq.'frgtyp'.or.trim(keywd(l4,k)(2:7)).eq.'FRGTYP')key_frag=1
if(trim(keywd(l4,k)(2:7)).eq.'frgtyp')key_frag=1
!if(trim(keywd(l4,k)(2:4)).eq.'nae'.or.trim(keywd(l4,k)(2:4)).eq.'NAE')then
if(trim(keywd(l4,k)(2:4)).eq.'nae')then
ll=len(trim(keywd(l4,k)))
do k1=1,ll+1
line5=keywd(l4,k)(k1:k1)
if(line5.eq.'=')then
line3=keywd(l4,k)(k1+1:ll)
read(line3,'(I2)')nae
endif
enddo
endif
enddo



do k=1,l3
!if(trim(keywd(l4,k)(2:4)).eq.'nao'.or.trim(keywd(l4,k)(2:4)).eq.'NAO')then
if(trim(keywd(l4,k)(2:4)).eq.'nao')then
ll=len(trim(keywd(l4,k)))
do k1=1,ll+1
line5=keywd(l4,k)(k1:k1)
if(line5.eq.'=')then
line3=keywd(l4,k)(k1+1:ll)
read(line3,'(I2)')nao
endif
enddo
endif
enddo


do k=1,l3
!if(trim(keywd(l4,k)(2:5)).eq.'nmul'.or.trim(keywd(l4,k)(2:5)).eq.'NMUL')then
if(trim(keywd(l4,k)(2:5)).eq.'nmul')then
ll=len(trim(keywd(l4,k)))
do k1=1,ll+1
line5=keywd(l4,k)(k1:k1)
if(line5.eq.'=')then
line3=keywd(l4,k)(k1+1:ll)
read(line3,'(I2)')mult
endif
enddo
endif
enddo

do k=1,l3
!if(trim(keywd(l4,k)(2:7)).eq.'chinst'.or.trim(keywd(l4,k)(2:7)).eq.'CHINST')then
if(trim(keywd(l4,k)(2:7)).eq.'chinst')then
ll=len(trim(keywd(l4,k)))
do k1=1,ll+1
line5=keywd(l4,k)(k1:k1)
if(line5.eq.'=')then
line3=keywd(l4,k)(k1+1:ll)
read(line3,'(I2)')chinst
endif
enddo
if(chinst.eq.1)then
flg1=0
itb=1
syb=3
nnb=2
radical=0
mnbond=0
endif
endif
enddo


do k=1,l3
!if(trim(keywd(l4,k)(2:4)).eq.'str'.or.trim(keywd(l4,k)(2:4)).eq.'STR')then
if(trim(keywd(l4,k)(2:4)).eq.'str')then
ll=len(trim(keywd(l4,k)))
do k1=1,ll+1
line5=keywd(l4,k)(k1:k1)
if(line5.eq.'=')then
line3=keywd(l4,k)(k1+1:ll)
!if(trim(line3).eq.'full'.or.trim(line3).eq.'FULL')flgst=1
if(trim(line3).eq.'full')flgst=1
!if(trim(line3).eq.'cov'.or.trim(line3).eq.'COV')flgst=2
if(trim(line3).eq.'cov')flgst=2
!if(trim(line3).eq.'ion'.or.trim(line3).eq.'ION')flgst=3
if(trim(line3).eq.'ion')flgst=3

endif
enddo
endif
enddo
rewind(21)

if(sl1(i).eq.iorb)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! reading of orbital section !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!print*,'orbital section start'

ll=0
l=0
do j=sl(i)+1,sl(i+1)-1
ll=len(trim(charst(j)))
line3=trim(charst(j)(1:1))
if(ll.eq.0)goto 300
if(trim(line3).eq.'#')goto 300
l=l+1
if(l.eq.1)goto 300
k2=0
do k=1,ll+1
line3=trim(charst(j)(k:k))
if(line3.eq.';')goto 301
if(line3.eq.'')then
k2=k2+1
gap2(k2)=k
endif
enddo
301 if(gap2(1).eq.1)then
k5=0
do k=2,k2
if(abs(gap2(k)-gap2(k-1)).ne.1)then
line3=trim(charst(j)(gap2(k-1)+1:gap2(k)-1))
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
if(line2.eq.'-')then
line5=trim(charst(j)(gap2(k-1)+1:gap2(k-1)+k3-1))
line6=trim(charst(j)(gap2(k-1)+k3+1:gap2(k)-1))
read(line5,'(I10)')orbn1
read(line6,'(I10)')orbn2
do k4=orbn1,orbn2
k5=k5+1
orbs(l,k5)=k4
enddo
goto 444
endif
enddo
k5=k5+1
read(line3,'(I10)')orbs(l,k5)
endif
444 enddo
endif
if(gap2(k).ne.1)then
k5=0
do k=1,k2
if(k.eq.1)then
line3=trim(charst(j)(1:gap2(k)-1))
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
if(line2.eq.'-')then
line5=trim(charst(j)(1:k3-1))
line6=trim(charst(j)(k3+1:gap2(k)-1))
read(line5,'(I10)')orbn1
read(line6,'(I10)')orbn2
do k4=orbn1,orbn2
k5=k5+1
orbs(l,k5)=k4
enddo
goto 445
endif
 enddo
k5=k5+1
read(line3,'(I10)')orbs(l,k5)
endif
if(k.ne.1)then
if(abs(gap2(k)-gap2(k-1)).ne.1)then
line3=trim(charst(j)(gap2(k-1)+1:gap2(k)-1))
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
if(line2.eq.'-')then
line5=trim(charst(j)(gap2(k-1)+1:gap2(k-1)+k3-1))
line6=trim(charst(j)(gap2(k-1)+k3+1:gap2(k)-1))
read(line5,'(I10)')orbn1
read(line6,'(I10)')orbn2
do k4=orbn1,orbn2
k5=k5+1
orbs(l,k5)=k4
enddo
goto 445
endif
 enddo
k5=k5+1
read(line3,'(I10)')orbs(l,k5)
endif
endif
445 enddo
endif

orbn(l)=k5
300 enddo
rewind(21)
orbsl=l
niao=l-1-nao
endif

do l=1,orbsl
!print*,'orbitals',orbsl,niao,(orbs(l,k),k=1,orbn(l))
print*,'orbitals',orbsl,niao,orbn(l),(orbs(l,k),k=1,4)
enddo




!do l=1,nfrag
!print*,'fragments',nfrag,(frag_at(l,k),k=1,frag_atn(l))
!enddo



if(sl1(i).eq.4)then

ll=0
l3=0
do j=sl(i)+1,sl(i+1)-1
l3=l3+1
gap1(1)=0
ll=len(trim(charst(j)))
!print*,'trim(charst(j))',trim(charst(j))
l=1
l1=0
do k=1,ll+1
line3=trim(charst(j)(k:k))
if(line3.eq.''.or.line3.eq.'-')then
l=l+1
gap1(l)=k
if(line3.eq.'-')dash(l)=l
endif
enddo
do k=2,l
!print*,'trim(charst(j)(gap1(k-1):gap1(k)))',gap1(k)-gap1(k-1),trim(charst(j)(gap1(k-1)+1:gap1(k)-1))
if(gap1(k)-gap1(k-1).gt.1)then
l1=l1+1
line3=trim(charst(j)(gap1(k-1)+1:gap1(k)-1))
!print*,'l3,l1',l3,l1
read(line3,'(I3)')bfimat(l3,l1)
!print*,'bfimat',bfimat(l3,l1)
if(dash(k).eq.k)then
l1=l1+1
bfimat(l3,l1)=0
endif
endif
enddo
bfiln(l3)=l1
enddo
!do k=1,l3
!print*,'bfimat(k,l2)',l1,(bfimat(k,k1),k1=1,bfiln(k))
!enddo
lastlbfi=l3
l2=0
j=lastlbfi
!print*,'l1l1l1',l1,j,bfiln(j)
do k=1,bfiln(j)-1
!print*,'kkkkkk',k,bfimat(j,k)
if(k.eq.1.and.bfimat(j,k).gt.1)then
do k1=1,bfimat(j,k)-1
l2=l2+1
freezorb(l2)=k1
!print*,'freezorb(l2)**',l2,freezorb(l2)
enddo
endif
if(bfimat(j,k+1).ne.0.or.k.ne.1)then
if(abs(bfimat(j,k+1)-bfimat(j,k)).ne.1)then
if(abs(bfimat(j,k+1)-bfimat(j,k)).ne.bfimat(j,k+1))then
!print*,'bfimat(j,k)+1,bfimat(j,k+1)',bfimat(j,k),bfimat(j,k+1),abs(bfimat(j,k+1)-bfimat(j,k))
do k1=bfimat(j,k)+1,bfimat(j,k+1)-1
l2=l2+1
freezorb(l2)=k1
!print*,'freezorb(l2)**',l2,freezorb(l2)
enddo
endif
endif
endif
enddo

frzn=l2

!do k1=1,frzn
!print*,'freezorb(l2)',frzn,freezorb(k1)
!enddo

endif
if(sl1(i).eq.ifrag)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! reading of fragment section !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!print*,'fragment section start'
k6=0
k7=0
k8=0
k9=0
do j=1,20
gap1(j)=0
enddo
ll=0
l=0
!print*,'sl',sl(i),sl(i+1)-1
do j=sl(i)+1,sl(i+1)-1
ll=len(trim(charst(j)))
line3=trim(charst(j)(1:1))
if(ll.eq.0)goto 310
if(trim(line3).eq.'#')goto 310
l=l+1
if(l.eq.1)goto 310
k2=0

do k=1,ll+1
line3=trim(charst(j)(k:k))
if(line3.eq.';')goto 311
if(line3.eq.'')then
k2=k2+1
gap2(k2)=k
endif
enddo
311 if(gap2(1).eq.1)then

if(key_frag.eq.1)then
do k=2,k2

if(abs(gap2(k)-gap2(k-1)).ne.1)then
line3=trim(charst(j)(gap2(k-1)+1:gap2(k)-1))

!print*,'line3:xmi',line3

do l1=orbsl-nao+1,orbsl
if(orbs(l1,1).eq.l-1)then
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
!if(line2.eq.'S'.or.line2.eq.'s')then
if(line2.eq.'s')then
k6=k6+1
orbsym(4,k6)=l1-1
goto 208
endif
enddo
do k3=1,ll
line2=trim(line3(k3:k3+1))
!if(line2.eq.'PX'.or.line2.eq.'Px'.or.line2.eq.'px')then
if(line2.eq.'px')then
k7=k7+1
orbsym(1,k7)=l1-1
goto 208
endif
enddo
do k3=1,ll
line2=trim(line3(k3:k3+1))
!if(line2.eq.'PY'.or.line2.eq.'Py'.or.line2.eq.'py')then
if(line2.eq.'py')then
k8=k8+1
orbsym(2,k8)=l1-1
goto 208
endif
enddo
do k3=1,ll
line2=trim(line3(k3:k3+1))
!if(line2.eq.'PZ'.or.line2.eq.'Pz'.or.line2.eq.'pz')then
if(line2.eq.'pz')then
k9=k9+1
orbsym(3,k9)=l1-1
goto 208
endif
enddo
endif
208 enddo
do l1=orbsl-nao+1,orbsl
if(orbs(l1,1).eq.l-1)then
l3=0
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
!if(line2.eq.'S'.or.line2.eq.'s')then
if(line2.eq.'s')then
l3=l3+1
norbsym1(l1-1,l3)='S'
!print*,'norbsym1(l1-1,l3)S',l1-1,norbsym1(l1-1,l3)
!print*,'k6k6',k6,orbsym(4,k6)
goto 218
endif
enddo
218 do k3=1,ll
line2=trim(line3(k3:k3+1))
!print*,'line2**',line2
!if(line2.eq.'PX'.or.line2.eq.'Px'.or.line2.eq.'px')then
if(line2.eq.'px')then
l3=l3+1
norbsym1(l1-1,l3)='X'
!print*,'norbsym1(l1-1,l3)X',l1-1,norbsym1(l1-1,l3)
!print*,'k7k7',k7
goto 219
endif
enddo
219 do k3=1,ll
line2=trim(line3(k3:k3+1))
!if(line2.eq.'PY'.or.line2.eq.'Py'.or.line2.eq.'py')then
if(line2.eq.'py')then
l3=l3+1
norbsym1(l1-1,l3)='Y'
!print*,'norbsym1(l1-1,l3)Y',l1-1,norbsym1(l1-1,l3)
!print*,'k8k8',k8
goto 220
endif
enddo
220 do k3=1,ll
line2=trim(line3(k3:k3+1))
!if(line2.eq.'PZ'.or.line2.eq.'Pz'.or.line2.eq.'pz')then
if(line2.eq.'pz')then
!print*,'line2inside',line2
l3=l3+1
norbsym1(l1-1,l3)='Z'
!print*,'norbsym1(l1-1,l3)Z',l1-1,norbsym1(l1-1,l3)
!print*,'k9k9',k9,orbsym(3,k9)
goto 221
endif
enddo
endif
221 num_norbsym1(l1-1)=l3
enddo
goto 210
endif
enddo
endif

210 k5=0
do k=2,k2-1
if(abs(gap2(k)-gap2(k-1)).ne.1)then
line3=trim(charst(j)(gap2(k)+1:gap2(k+1)-1))
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
if(line2.eq.'-')then
line5=trim(charst(j)(gap2(k)+1:gap2(k)+k3-1))
line6=trim(charst(j)(gap2(k)+k3+1:gap2(k+1)-1))
read(line5,'(I10)')orbn1
read(line6,'(I10)')orbn2
do k4=orbn1,orbn2
k5=k5+1
frag_at(l,k5)=k4
!print*,'frag_at',k5,frag_at(l,k5)
enddo
goto 454
endif
enddo
!print*,'k55555',k5
k5=k5+1
read(line3,'(I10)')frag_at(l,k5)
endif
454 enddo
!print*,'k5555',k5
frag_atn(l)=k5
!print*,'sourav_frag_atn',frag_atn(l)
endif
if(gap2(1).ne.1)then
if(key_frag.eq.1)then
line3=trim(charst(j)(1:gap2(1)-1))
!print*,'line3:xmi',line3

do l1=orbsl-nao+1,orbsl
if(orbs(l1,1).eq.l-1)then
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
!print*,'line2**',line2,line3
!if(line2.eq.'S'.or.line2.eq.'s')then
if(line2.eq.'s')then
k6=k6+1
orbsym(4,k6)=l1-1
!print*,'k6k6',k6,orbsym(4,k6)
goto 209
endif
enddo
do k3=1,ll
line2=trim(line3(k3:k3+1))
!print*,'line2**',line2
!if(line2.eq.'PX'.or.line2.eq.'Px'.or.line2.eq.'px')then
if(line2.eq.'px')then
k7=k7+1
orbsym(1,k7)=l1-1
!print*,'k7k7',k7
goto 209
endif
enddo
do k3=1,ll
line2=trim(line3(k3:k3+1))
!if(line2.eq.'PY'.or.line2.eq.'Py'.or.line2.eq.'py')then
if(line2.eq.'py')then
k8=k8+1
orbsym(2,k8)=l1-1
!print*,'k8k8',k8
goto 209
endif
enddo
do k3=1,ll
line2=trim(line3(k3:k3+1))
!if(line2.eq.'PZ'.or.line2.eq.'Pz'.or.line2.eq.'pz')then
if(line2.eq.'pz')then
!print*,'line2inside',line2
k9=k9+1
orbsym(3,k9)=l1-1
!print*,'k9k9',k9,orbsym(3,k9)
goto 209
endif
enddo
endif
!enddo
209 enddo
do l1=orbsl-nao+1,orbsl
if(orbs(l1,1).eq.l-1)then
l3=0
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
!print*,'line2**',line2
!if(line2.eq.'S'.or.line2.eq.'s')then
if(line2.eq.'s')then
l3=l3+1
norbsym1(l1-1,l3)='S'
!print*,'norbsym1(l1-1,l3)S',l1-1,norbsym1(l1-1,l3)
!print*,'k6k6',k6,orbsym(4,k6)
goto 318
endif
enddo
318 do k3=1,ll
line2=trim(line3(k3:k3+1))
!print*,'line2**',line2
!if(line2.eq.'PX'.or.line2.eq.'Px'.or.line2.eq.'px')then
if(line2.eq.'px')then
l3=l3+1
norbsym1(l1-1,l3)='X'
print*,'norbsym1(l1-1,l3)X',l1-1,norbsym1(l1-1,l3)
!print*,'k7k7',k7
goto 319
endif
enddo
319 do k3=1,ll
line2=trim(line3(k3:k3+1))
!if(line2.eq.'PY'.or.line2.eq.'Py'.or.line2.eq.'py')then
if(line2.eq.'py')then
l3=l3+1
norbsym1(l1-1,l3)='Y'
!print*,'norbsym1(l1-1,l3)Y',l1-1,norbsym1(l1-1,l3)
!print*,'k8k8',k8
goto 320
endif
enddo
320 do k3=1,ll
line2=trim(line3(k3:k3+1))
!if(line2.eq.'PZ'.or.line2.eq.'Pz'.or.line2.eq.'pz')then
if(line2.eq.'pz')then
!print*,'line2inside',line2
l3=l3+1
norbsym1(l1-1,l3)='Z'
!print*,'norbsym1(l1-1,l3)Z',l1-1,norbsym1(l1-1,l3)
!print*,'k9k9',k9,orbsym(3,k9)
goto 321
endif
enddo
endif
321 num_norbsym1(l1-1)=l3
enddo
endif

k5=0
do k=2,k2
!print*,'k2222',k,k2,charst(j)
line3=trim(charst(j)(gap2(k-1)+1:gap2(k)-1))
if(trim(line3).eq.'') goto 404
!print*,'fragk=2',line3
ll=len(trim(line3))
do k3=1,ll
line2=trim(line3(k3:k3))
if(line2.eq.'-')then
line5=trim(charst(j)(gap2(k-1)+1:gap2(k-1)+k3-1))
line6=trim(charst(j)(gap2(k-1)+k3+1:gap2(k)-1))
read(line5,'(I10)')orbn1
read(line6,'(I10)')orbn2
do k4=orbn1,orbn2
k5=k5+1
frag_at(l,k5)=k4
enddo
goto 455
endif
 enddo
k5=k5+1
404 read(line3,'(I10)')frag_at(l,k5)
455 enddo
frag_atn(l)=k5
!print*,'sourav_frag_atnnnn',frag_atn(l)
endif
310 enddo
rewind(21)
nfrag=l

!do i1=1,nfrag
!print*,'frag_atoms',nfrag,(frag_at(i1,j),j=1,frag_atn(i1))
!enddo
norbsym(4)=k6
norbsym(1)=k7
norbsym(2)=k8
norbsym(3)=k9


endif

enddo


vacorb=nao-nae
nlp=nae-nao
if(vacorb.gt.1)nlp=0
nlast=(mult-1)

if(chinst.eq.1)call input
call read_info

!print*,'exit read_xmi'
500 return
end subroutine read_xmi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

integer::i,i1,i2,i3,k,k1,k2,k3,l,j,unkwn&
,nbasatom1(20),cnt(20),ufrzorbn,ovlp_basis_num(20,10),MNP,iO
integer::fragorb(30),frag_cent(50),frag_cent1(50),k4,k5,nfragcent,n,k9,k10,linenumi,linenumj
integer::bfiexists,n1,n2,n3,n4,l1,l2,j1,i4,sig(20),l3,l4,act_at_num(20),val_state_num(20) &
,pi1(20),pi2(20),k6,k7,nsig,npi1,npi2,k8,ovlp_at_bas_num(20,10),nbas_slab(100)
real*8::oval,a(4)
!biasval(20),least
character(len=5)::atbas(500),act_orb_typ(20,10),atbas_frz(50),line
logical :: fileexists
integer::atsymset(20,20),nsym,syn(50),at_sym(50)
real*8::coordx(100),coordy(100),coordz(100)

common/ats/atsymset,nsym,syn,at_sym
common/coordinate/coordx,coordy,coordz

print*,'enter read_info'

ovlp_int=0
INQUIRE(FILE='x1e.int',EXIST=fileexists)
IF (fileexists) THEN
open (unit=29,file='x1e.int',status='old')
read(29,*)ufrzorbn
ovlp_int=0
!close(29)

!open (unit=29,file='x1e.int',status='old')
MNP=0
do
read(29,'(a)',iostat=io)

if(io.ne.0)exit
MNP=MNP+1
enddo
rewind(29)


do i=1,MNP
read(29,*),line
if(trim(line(1:5)).eq.'SSF0')lnum=i
enddo

!close(29)
ELSE
PRINT*,' x1e.int file does not exist as a result I could not able &
to check if your system has frozen orbitals or not hope you have previded &
correct information otherwise result will be wrong'
ufrzorbn=0
!stop
ENDIF


INQUIRE(FILE='INFO',EXIST=fileexists)
IF (fileexists) THEN
open (unit=39,file='INFO',status='old')
ELSE
PRINT*,'SORRY INFO file does not exist'
stop
ENDIF
read(39,*),tot_atom,nbasis,unkwn
k=0
do i=1,tot_atom
read(39,*),nbasatom(i)
k=k+nbasatom(i)
nbas_slab(i)=k
!print*,'nbasatom',nbasatom(i),nbas_slab(i)
enddo



bfiexists=0
do i=1,10
if(sl1(i).eq.4)bfiexists=1
enddo

if(ufrzorbn.ne.0)then
if(ufrzorbn.lt.k.and.bfiexists.eq.0)then
write(*,*)'it seems that you have frozen orbitals but you forget to put the $bfi section'
stop
else

k=0
do i=1,tot_atom
l=0
k=k+nbasatom(i)
cnt(i)=k
do j=1,frzn
if(freezorb(j).le.cnt(i))then
if(i.gt.1)then
if(freezorb(j).le.cnt(i-1))goto 100 
endif
l=l+1
endif
100 enddo
nbasatom(i)=nbasatom(i)-l
enddo

endif
endif

if(key_frag.eq.0)then
do i=1,200
do j=1,20
atoset(i,j)=0
enddo
enddo

!do i=1,orbsl
!print*,'orbn**********',tot_atom,orbn(i),(orbs(i,j),j=1,orbn(i))
!enddo

!do i=2,niao+1
!j=0
!k=0
!!print*,'iiiiii',i,(orbs(i,l),l=1,orbn(i))
!do i1=1,tot_atom
!!print*,'nbas',nbasatom(i1)
!j=j+nbasatom(i1)
!do l=1,orbn(i)
!if(orbs(i,l).gt.j.or.orbs(i,l).le.j-nbasatom(i1))goto 310
!k=k+1
!atst_inact(i-1,k)=i1
!goto 311
!310 enddo
!311 enddo
!atst_inact_nm(i-1)=k
!atst_inact_val(i-1)=2/k
!enddo
!!k=0
!!do i=1,niao
!!print*,(atst_inact(i,k),k=1,atst_inact_nm(i))
!!enddo
!do i1=1,tot_atom
!inact_elecs(i1)=0
!enddo
!do i1=1,tot_atom
!do i=1,niao
!do i2=1,atst_inact_nm(i)
!if(i1.eq.atst_inact(i,i2))then
!inact_elecs(i1)=inact_elecs(i1)+atst_inact_val(i)
!endif
!enddo
!enddo
!enddo

!print*,(inact_elecs(i),i=1,tot_atom)
!stop

k5=0
j=0
k1=0
k3=0
k4=0
do i1=1,tot_atom
k=0
j=j+nbasatom(i1)
do i=orbsl-nao+1,orbsl
!print*,'active_orbs',i,orbn(i),orbsl,nao,mult
k2=0
do l=1,orbn(i)
if(orbs(i,l).gt.j.or.orbs(i,l).le.j-nbasatom(i1))goto 200
k2=k2+1
enddo
k5=k5+1
k=k+1

do i4=1,k3
if(i1.eq.active_atoms(i4))goto 207
enddo
k3=k3+1
active_atoms(k3)=i1
print*,'active_atoms_1',active_atoms(k3)

207 atoset(i1,k)=i-1
atn(i1)=k 
active_orbs(k5)=i-1
atm_nb_orbs(k5)=i1
!print*,'orbs(i),j**',orbs(i,l),j,i,i1,k
goto 201
200 if(k2.ne.0.and.num_frag_cntr.eq.0)then
k1=k1+1
fragorb(k1)=i-1
if(orbs(i,1).le.j.and.orbs(i,1).gt.j-nbasatom(i1))then
k5=k5+1
k=k+1
k4=k4+1

!print*,'k3k3k3**',k3
do i4=1,k3
if(i1.eq.active_atoms(i4))goto 206
!print*,'i4i4i4',i4
enddo
k3=k3+1
!print*,'k3k3**',k3
active_atoms(k3)=i1
print*,'active_atoms_2',active_atoms(k3)

206 frag_cent(k4)=i1
atoset(i1,k)=i-1
atn(i1)=k 
active_orbs(k5)=i-1
atm_nb_orbs(k5)=i1
endif
endif
if(k2.ne.0.and.num_frag_cntr.ne.0)then
do i3=1,num_frag_cntr
j1=0
do i2=1,frag_cntr(i3)
j1=j1+nbasatom(i2)
enddo
!print*,'j1',j1,frag_cntr(i3)
!if(orbs(i,1).le.j1.and.orbs(i,1).gt.j1-nbasatom(i2))then
k5=k5+1
k=k+1
k4=k4+1
!print*,'k3k3k3:frag_cntr',k3
do i4=1,k3
print*,'i3i3i3',i3,active_atoms(i4)
if(frag_cntr(i3).eq.active_atoms(i4))goto 205
!print*,'i4i4i4',i4
enddo
k3=k3+1
!print*,'k3k3sourav',k3
active_atoms(k3)=frag_cntr(i3)
print*,'active_atoms_3',active_atoms(k3)

205 frag_cent(k4)=frag_cntr(i3)
atoset(frag_cntr(i3),k)=i-1
active_orbs(k5)=i-1
atm_nb_orbs(k5)=frag_cntr(i3)
atn(frag_cntr(i3))=k
!print*,'kkkatn',k
!endif

enddo
endif
201 enddo
!print*,'k33333',k3,i1
enddo

atom=k3
nactorb=k5
!print*,'active_orbs',nactorb,(active_orbs(i),i=1,nactorb)
!print*,'atn',atom,(atn(active_atoms(i)),i=1,atom)
print*,'active_atom',atom,(active_atoms(i),i=1,atom)
!print*,'atm_nb_orbs(k5)',(atm_nb_orbs(i),i=1,nactorb)
do i=1,tot_atom
print*,'atoset:info',(atoset(i,j),j=1,atn(i))
!
enddo

k5=1
frag_cent1(k5)=frag_cent(1)
do i=2,k4
do j=1,k5
if(frag_cent(i).eq.frag_cent1(j))goto 320
enddo
k5=k5+1
frag_cent1(k5)=frag_cent(i-1)
320 enddo
nfragcent=k5
!print*,'num of fragorbs',(fragorb(k),k=1,k1)
!print*,'frag_cent1',(frag_cent1(i),i=1,k5)

if(k1.ne.0.and.num_frag_cntr.eq.0)then
if(k1.eq.1)print*,'You have fragment orbital ',fragorb(1),' . Please mention the &  
central atoms of the fragmants or start the orbitals with basis function of the central atom'
if(k1.gt.1)print*,'You have fragment orbitals ',(fragorb(k),k=1,k1),' . Please mention the &  
central atoms of the fragmants or start the orbitals with basis function of the central atom'
endif

l2=0
do k1=1,tot_atom
l1=0
if(atoset(k1,1).eq.0)goto 321
l2=l2+1
do k2=1,20
if(atoset(k1,k2).ne.0)then
l1=l1+1
!atoset(l2,l1)=atommat(k1,k2)
!print*,'atoset',(atoset(k1,k2),k2=1,5)
endif
enddo
atn(k1)=l1
321 enddo
atom=l2

!do i=1,tot_atom
!!print*,'iii',i
!print*,'atoset_info_2',(atoset(i,j),j=1,5)
!enddo
do i=1,tot_atom
!!print*,'iii**',i
print*,'atoset_info_2',(atoset(i,j),j=1,atn(i))
enddo

endif

do i=1,tot_atom
read(39,*),symat(i),symatno(i),coordx(i),coordy(i),coordz(i)
coordx(i)=coordx(i)*0.529177
coordy(i)=coordy(i)*0.529177
coordz(i)=coordz(i)*0.529177
!print*,'coord**',symat(i),symatno(i),coordx(i),coordy(i),coordz(i)
all_at_num(i)=int(symatno(i))
!print*,'all_at_num',all_at_num(i)
enddo
if(key_frag.eq.0)then
do i=1,atom
act_at_num(i)=symatno(active_atoms(i))
val_state_num(i)=valence_state(act_at_num(i))
enddo

!print*,'act_at_num',(act_at_num(i),i=1,atom)
!print*,'val_state_num',(val_state_num(i),i=1,atom)
endif


!if(k1.ne.0)then
!write(*,*)'as you did not put the central atom of the fragments we took the atom number ='&
!,(frag_cent1(i),i=1,k5),'(',(symat(frag_cent1(i)),i=1,nfragcent),')','as fragment center'
!endif
!990 format(a,1x,I2,1x,a,a,a,a)
!print*,'tot_atom:info',tot_atom
!n=0
!n1=0
!do i=1,tot_atom
!if(atoset(i,1).eq.0)goto 189
!print*,'**********',atoset(i,1),i,symat(i)
!n1=n1+1
!do i3=1,88
!if(symat(i).eq.at_list(i3))then
!print*,'symat(i).eq.at_list(i3)',symat(i),at_list(i3)
!n=n+1
!actv_atom(n)=i
!at_num(n)=i3
!endif
!enddo
!189 enddo
!if(n1.ne.n)then
!n=0
!do i=1,tot_atom
!if(atoset(i,1).eq.0)goto 190
!print*,'**********',atoset(i,1),i,symat(i)
!do i3=1,88
!if(symat(i).eq.at_list_bold(i3))then
!print*,'symat(i).eq.at_list(i3)',symat(i),at_list_bold(i3)
!n=n+1
!actv_atom(n)=i
!at_num(n)=i3
!endif
!enddo
!190 enddo
!
!endif
!do i=1,n
!print*,'at_num',symat(i),at_num(i),actv_atom(i)
!enddo

n=0
do i=1,nbasis
if(frzn.ne.0)then
do j=1,frzn
!print*,'frzn',frzn,freezorb(j)
if(freezorb(j).eq.i)then
read(39,*)atbas_frz(i)
all_atbas(i)=atbas_frz(i)
goto 197
endif
enddo
endif
n=n+1
read(39,*),atbas(n)
all_atbas(i)=atbas(n)
!print*,'atbas:sourav',i,n,atbas(n)
197 enddo

!do i=1,nbasis
!print*,'all_atbas(i)',i,all_atbas(i)
!enddo

!enddo


if(key_frag.eq.0)then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! determination of the type of the orbitals 'pi or sigma' !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
norbsym(4)=0
norbsym(1)=0
norbsym(2)=0
norbsym(3)=0

k=0
k5=0
k8=0
k4=0
k6=0
k1=0
k3=0
do i=orbsl-nao+1,orbsl
k7=0
n=0
!print*,'atom',atom
do i2=1,atom
do i3=1,atn(active_atoms(i2))
print*,'ovlp_act_atm,atoset',active_atoms(i2),atoset(active_atoms(i2),i3)
if(atoset(active_atoms(i2),i3).eq.i-1)then
print*,'active_atoms(i2)',active_atoms(i2),nbasatom(3)
do j=1,active_atoms(i2)
k7=k7+1
n=n+nbasatom(j)
enddo
endif
enddo
enddo
!print*,'n,nbasatom(k7)',k7,n,nbasatom(k7),n-nbasatom(k7),i-1,orbn(i)
!do i2=1,orbn(i)
!print*,'test_ovlp_atbas(orbs(i,i2)',i,orbs(i,i2),atbas(orbs(i,i2))
!enddo




!do k10=1,num_act_orb_typ(i)
!l4=0
!do k9=1,nbasis
!if(all_atbas(k9).eq.act_orb_typ(i,k10))then
!l4=l4+1
!if(l4.eq.valence_state((i)))then
!ovlp_basis_num(i,k10)=k9
!endif
!endif
!enddo
!enddo



do i2=1,orbn(i)
if(orbs(i,i2).le.n.and.orbs(i,i2).gt.n-nbasatom(k7))then
if(atbas(orbs(i,i2)).eq.'S')then
k=k+1
orbsym(4,k)=i-1
!print*,'orbsym(4,k)',orbsym(4,k)
goto 208
endif
endif
enddo
do i2=1,orbn(i)
if(orbs(i,i2).le.n.and.orbs(i,i2).gt.n-nbasatom(k7))then
if(atbas(orbs(i,i2)).eq.'X')then
k3=k3+1
orbsym(1,k3)=i-1
!print*,'orbsym(1,k)',orbsym(1,k3)
goto 208
endif
endif
enddo
do i2=1,orbn(i)
if(orbs(i,i2).le.n.and.orbs(i,i2).gt.n-nbasatom(k7))then
if(atbas(orbs(i,i2)).eq.'Y')then
k4=k4+1
orbsym(2,k4)=i-1
!print*,'orbsym(2,k)',orbsym(2,k4)
goto 208
endif
endif
enddo
do i2=1,orbn(i)
if(orbs(i,i2).le.n.and.orbs(i,i2).gt.n-nbasatom(k7))then
!print*,'atbas(orbs(i,i2))',atbas(orbs(i,i2))
if(atbas(orbs(i,i2)).eq.'Z')then
k8=k8+1
orbsym(3,k8)=i-1
!print*,'orbsym(3,k)',orbsym(3,k8)
goto 208
endif
endif
enddo

208 enddo
norbsym(4)=k
norbsym(1)=k3
norbsym(2)=k4
norbsym(3)=k8


if(ovlp_int.eq.1.and.ovopt.eq.1)then
n2=0
do i=1,atom
n=0
n1=0
n4=0
do i2=1,active_atoms(i)-1
n=n+nbasatom(i2)
enddo
do i2=1,active_atoms(i)
n1=n1+nbasatom(i2)
enddo
do i2=1,frzn
if(freezorb(i2).ge.n.and.freezorb(i2).le.n1)then
n4=n4+1
endif
enddo
n2=n2+n4
do i3=1,atn(active_atoms(i))
l1=0
l2=0
l3=0
l4=0
k9=0
if(n4.eq.0)then
do i2=1,orbn(atoset(active_atoms(i),i3)+1)
print*,'orbn(i)',i2,orbn(atoset(active_atoms(i),i3)+1),orbs(atoset(active_atoms(i),i3)+1,i2)+n2
n3=orbs(atoset(active_atoms(i),i3)+1,i2)+n2
if(n3.lt.n.or.n3.gt.n1)goto 531
if(all_atbas(n3).eq.'S')then
l1=l1+1
if(l1.eq.val_state_num(i))goto 530
endif
if(all_atbas(n3).eq.'X')then
l2=l2+1
if(l2.eq.val_state_num(i))goto 530
endif
if(all_atbas(n3).eq.'Y')then
l3=l3+1
if(l3.eq.val_state_num(i)-1)goto 530
endif
if(all_atbas(n3).eq.'Z')then
l4=l4+1
if(l4.eq.val_state_num(i)-1)goto 530
endif
goto 531
530 k9=k9+1
act_orb_typ(atoset(active_atoms(i),i3),k9)=all_atbas(n3)
ovlp_bas_num(atoset(active_atoms(i),i3),k9)=n3
531 enddo
endif
if(n4.ne.0)then
do i2=1,orbn(atoset(active_atoms(i),i3)+1)
!print*,'orbn(i)',i2,orbn(atoset(active_atoms(i),i3)+1),orbs(atoset(active_atoms(i),i3)+1,i2)+n2
n3=orbs(atoset(active_atoms(i),i3)+1,i2)+n2
!print*,'all_atbas(n3)',n3,all_atbas(n3)
if(n3.lt.n.or.n3.gt.n1)goto 631
if(all_atbas(n3).eq.'S'.and.l1.eq.0)then
l1=1
goto 630
endif
if(all_atbas(n3).eq.'X'.and.l2.eq.0)then
l2=1
goto 630
endif
if(all_atbas(n3).eq.'Y'.and.l3.eq.0)then
l3=1
goto 630
endif
if(all_atbas(n3).eq.'Z'.and.l4.eq.0)then
l4=1
goto 630
endif
goto 631
630 k9=k9+1
act_orb_typ(atoset(active_atoms(i),i3),k9)=all_atbas(n3)
ovlp_bas_num(atoset(active_atoms(i),i3),k9)=n3
631 enddo
endif
num_act_orb_typ(atoset(active_atoms(i),i3))=k9

enddo
!do i3=1,atn(active_atoms(i))
!print*,'act_orb_typ',atoset(active_atoms(i),i3),(act_orb_typ(atoset(active_atoms(i),i3),k10),k10=1,k9)&
!,(ovlp_bas_num(atoset(active_atoms(i),i3),k10),k10=1,k9)
!enddo

enddo

endif

endif

print*,'norbsymmmmm',norbsym(1),norbsym(2),norbsym(3),norbsym(4),atom
!if(ovlp_int.eq.1.and.ovopt.eq.1)then
!do i1=1,nactorb-1
!do j=i1+1,nactorb
!!do i1=1,nactorb
!!do j=1,nactorb
!oval=0.0
!do i2=1,num_act_orb_typ(active_orbs(i1))
!do i3=1,num_act_orb_typ(active_orbs(j))
!n1=ovlp_bas_num(active_orbs(i1),i2)
!n2=ovlp_bas_num(active_orbs(j),i3)
!print*,'n1,n2,nbasis',n1,n2,nbasis,i2,i3
!linenumi=(nbasis*(n1-1)+n2)/4
!linenumj=mod(nbasis*(n1-1)+n2,4)
!n4=linenumj
!if(linenumj.eq.0)then
!linenumi=linenumi-1
!n4=4
!endif
!print*,'lnum',lnum,linenumi,linenumj
!rewind(29)
!do i=1,lnum
!read(29,*),line
!!print*,'line1',i,line
!enddo
!!if(linenumi.eq.0)goto 534
!do i=lnum+1,lnum+linenumi
!read(29,*)line
!!print*,'line2',i,line
!enddo
!534 read(29,*)(a(i),i=1,4)
!print*,(a(i),i=1,4)
!!read(29,*)line
!!print*,line
!n3=num_act_orb_typ(active_orbs(i1))*num_act_orb_typ(active_orbs(j))
!oval=oval+a(n4)/n3
!print*,'oval',oval,a(n4),n3,n4,a(n4)/n3
!rewind(29)
!enddo
!enddo
!orb_ovlp_mat1(i1,j)=oval
!enddo
!enddo
!
!do i1=1,nactorb
!orb_ovlp_mat1(i1,i1)=1.0
!enddo
!
!do i1=1,nactorb-1
!do j=i1+1,nactorb
!orb_ovlp_mat1(j,i1)=orb_ovlp_mat1(i1,j)
!enddo
!enddo
!
!do i1=1,nactorb
!print*,'orb_ovlp_mat1(i1,j)',(orb_ovlp_mat1(i1,j),j=1,nactorb)
!enddo
!
!endif

k5=0
k6=0
do i=1,4
if(norbsym(i).ne.0)then
k5=k5+1
if(i.lt.4)then
k6=1
else
k6=k6+1
endif
at_sym(k5)=k6
syn(k5)=norbsym(i)
do j=1,norbsym(i)
atsymset(k5,j)=orbsym(i,j)
enddo
endif
enddo
nsym=k5


do i=1,nsym
print*,'atsymset',nsym,at_sym(i),(atsymset(i,j),j=1,syn(i))
enddo
!print*,'sig',(sig(i),i=1,nsig)
!print*,'pi1',(pi1(i),i=1,npi1)
!print*,'pi2',(pi2(i),i=1,npi2)

!do i=1,nsym
!print*,'atsymset:info',at_sym(i),(atsymset(i,j),j=1,syn(i))
!enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print*,'exit read_info'
call geocal(coordx,coordy,coordz)
return
end subroutine read_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine geocal(coordx,coordy,coordz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

integer:: i,j,ij,j1,j2,j3,k,ii1,i1,i2,i3,i4,k1,k2,k3,k4,k9,k10,nland,nisland(50),islands(20,50),nncatm,ncatm,&
nnmat(100,100),nnact(100),nnact_1(100),active_orb2(20),active_orb1(20),l,l1,l2,l3,l4,val_st(20)&
,ind(100),act_at_num(20),val_state_num(100),n,n1,n2,n3,n4&
,val_state_num1(100),linenumi,linenumj,active_atoms_sort(100)
real*8::least,bond_dist,ddist,kk,island_mat(20,20)
real*8::coordx(100),coordy(100),coordz(100),oval,a(4)
character(len=5)::act_orb_typ(20,10),line

print*,'enter geocal'
sig_sym_flg=0
if(key_frag.eq.1)then
do i=1,200
do j=1,20
atoset(i,j)=0
enddo
enddo
k3=0
do i=1,tot_atom
ind(i)=1
enddo
do k=orbsl-nao+1,orbsl
if(frag_atn(orbs(k,1)).gt.1)then
write(*,990),'The orbital',k-1,'is a fragment orbital and we took the atom number'&
,frag_at(orbs(k,1)+1,1),'as the central atom.'
write(*,*),'If it is not right please put the central atom number at first place &
in the row in the $frag section and run again.'
endif
if(atoset(frag_at(orbs(k,1)+1,1),ind(frag_at(orbs(k,1)+1,1))).ne.0)&
ind(frag_at(orbs(k,1)+1,1))=ind(frag_at(orbs(k,1)+1,1))+1
if(num_frag_cntr.ne.0)then
do j=1,frag_atn(orbs(k,1)+1)
do i=1,num_frag_cntr
if(frag_cntr(i).eq.frag_at(orbs(k,1)+1,j))then
frag_at(orbs(k,1)+1,1)=frag_cntr(i)
endif
enddo
enddo
endif

atoset(frag_at(orbs(k,1)+1,1),ind(frag_at(orbs(k,1)+1,1)))=k-1
!print*,'atoset22',frag_at(orbs(k,1)+1,1),atoset(frag_at(orbs(k,1)+1,1),ind(frag_at(orbs(k,1)+1,1)))
!print*,'frag_at(orbs(k,1)+1,1)',frag_at(orbs(k,1)+1,1),k-1,k1,atoset(frag_at(orbs(k,1)+1,1),k1)
do i=1,8
print*,'atn:geo_1',atn(active_atoms(i)),active_atoms(i)
enddo
do i=1,k3
if(active_atoms(i).eq.frag_at(orbs(k,1)+1,1))goto 266
enddo
k3=k3+1
active_atoms(k3)=frag_at(orbs(k,1)+1,1)

266 enddo
do i=1,k3
atn(active_atoms(i))=ind(active_atoms(i))
act_at_num(i)=symatno(active_atoms(i))
val_state_num(i)=valence_state(act_at_num(i))
val_state_num1(i)=val_state_num(i)
enddo
do i=1,k3
print*,'atn:geo',atn(active_atoms(i)),active_atoms(i)
enddo
atom=k3

!if(ovlp_int.eq.1.and.ovopt.eq.1)then
!n2=0
!do i=1,atom
!n=0
!n1=0
!n4=0
!do i2=1,active_atoms(i)-1
!n=n+nbasatom(i2)
!enddo
!do i2=1,active_atoms(i)
!n1=n1+nbasatom(i2)
!enddo
!!print*,'frzn',frzn,freezorb(j)
!do i3=1,atn(active_atoms(i))
!l2=0
!l3=0
!l4=0
!k9=0
!!print*,'atn(active_atoms(i)):geo',i3,atn(active_atoms(i))
!do i2=1,num_norbsym1(atoset(active_atoms(i),i3))
!l1=0
!!print*,'orbn(i):geo',i2,orbn(atoset(active_atoms(i),i3)+1),orbs(atoset(active_atoms(i),i3)+1,i2)+n2
!do n3=1+n,n1
!!print*,'all_atbas(n3):geo',n3,all_atbas(n3),norbsym1(atoset(active_atoms(i),i3),i2)
!!if(n3.lt.n.or.n3.gt.n1)goto 531
!val_state_num1(i)=val_state_num(i)
!if(norbsym1(atoset(active_atoms(i),i3),i2).eq.'X')val_state_num1(i)=val_state_num1(i)-1
!if(norbsym1(atoset(active_atoms(i),i3),i2).eq.'Y')val_state_num1(i)=val_state_num1(i)-1
!if(norbsym1(atoset(active_atoms(i),i3),i2).eq.'Z')val_state_num1(i)=val_state_num1(i)-1
!!print*,'val_state_num(i)',val_state_num1(i)
!if(all_atbas(n3).eq.norbsym1(atoset(active_atoms(i),i3),i2))then
!l1=l1+1
!if(l1.eq.val_state_num1(i))goto 530
!endif
!enddo
!goto 531
!530 k9=k9+1
!act_orb_typ(atoset(active_atoms(i),i3),k9)=all_atbas(n3)
!ovlp_bas_num(atoset(active_atoms(i),i3),k9)=n3
!!print*,'ovlp_bas_num(atoset(active_atoms(i),i3),k9)i:geo',ovlp_bas_num(atoset(active_atoms(i),i3),k9),&
!!all_atbas(n3)
!531 enddo
!num_act_orb_typ(atoset(active_atoms(i),i3))=k9
!
!enddo
!!do i3=1,atn(active_atoms(i))
!!print*,'act_orb_typ:geo',atoset(active_atoms(i),i3),(act_orb_typ(atoset(active_atoms(i),i3),k10),k10=1,&
!!num_act_orb_typ(atoset(active_atoms(i),i3)))&
!!,(ovlp_bas_num(atoset(active_atoms(i),i3),k10),k10=1,k9)
!!enddo
!
!enddo
!
!endif

endif

l3=0
l2=0
do k1=1,tot_atom
l1=0
if(atoset(k1,1).eq.0)goto 326
l2=l2+1
do k2=1,20
if(atoset(k1,k2).ne.0)then
l1=l1+1
l3=l3+1
active_orbs(l3)=atoset(k1,k2)
atm_nb_orbs(l3)=k1
endif
enddo
atn(k1)=l1
326 enddo
atom=l2
nactorb=l3

!print*,'nactorb',nactorb,tot_atom
print*,'active_atom**',atom,(active_atoms(i),i=1,atom)
!print*,'atn:atn:sourav',(atn(active_atoms(i)),i=1,atom)
do i=1,nactorb
active_orb2(i)=active_orbs(i)
atm_nb_orbs1(i)=atm_nb_orbs(i)
enddo

k1=0
k=0
436 least=100
do i=1,nactorb
do j=1,k1
if(active_orb1(j).eq.i)goto 337
enddo
if(active_orbs(i).lt.least)then
least=active_orbs(i)
k=i
endif
337 enddo
k1=k1+1
active_orb1(k1)=k
if(k1.lt.nactorb)goto 436
do i=1,nactorb
active_orbs(i)=0
atm_nb_orbs(i)=0
enddo
do i=1,nactorb
active_orbs(i)=active_orb2(active_orb1(i))
atm_nb_orbs(i)=atm_nb_orbs1(active_orb1(i))
enddo



do i=1,tot_atom
print*,'atosetat_geo',atn(i),(atoset(i,j),j=1,atn(i))
enddo

do i=1,nactorb
print*,'sourav4',i
do j=1,4
do k=1,norbsym(j)
if(norbsym(j).ne.0)then
if(active_orbs(i).eq.orbsym(j,k))then
atm_nb_sym(i)=j
endif
endif
enddo
enddo
enddo
!open (unit=29,file='x1e.int',status='old')
!if(ovlp_int.eq.1.and.ovopt.eq.1)then
!print*,'i am here'
!stop
!do i1=1,nactorb-1
!do j=i1+1,nactorb
!oval=0.0
!do i2=1,num_act_orb_typ(active_orbs(i1))
!do i3=1,num_act_orb_typ(active_orbs(j))
!n1=ovlp_bas_num(active_orbs(i1),i2)
!n2=ovlp_bas_num(active_orbs(j),i3)
!linenumi=(nbasis*(n1-1)+n2)/4
!linenumj=mod(nbasis*(n1-1)+n2,4)
!n4=linenumj
!if(linenumj.eq.0)then
!linenumi=linenumi-1
!n4=4
!endif
!rewind(29)
!do i=1,lnum
!read(29,*),line
!enddo
!do i=lnum+1,lnum+linenumi
!read(29,*)line
!enddo
!534 read(29,*)(a(i),i=1,4)
!n3=num_act_orb_typ(active_orbs(i1))*num_act_orb_typ(active_orbs(j))
!oval=oval+a(n4)/n3
!rewind(29)
!enddo
!enddo
!orb_ovlp_mat1(i1,j)=oval
!enddo
!enddo
!
!print*,'sourav5'
!do i1=1,nactorb
!orb_ovlp_mat1(i1,i1)=1.0
!enddo
!
!do i1=1,nactorb-1
!do j=i1+1,nactorb
!orb_ovlp_mat1(j,i1)=orb_ovlp_mat1(i1,j)
!enddo
!enddo
!
!!do i1=1,nactorb
!!print*,'orb_ovlp_mat1(i1,j)',(orb_ovlp_mat1(i1,j),j=1,nactorb)
!!enddo
!endif
!close(29)

!do i=1,nactorb
!print*,'norbsym1(active_atom(i))',active_orbs(i),(norbsym1(active_orbs(i),i1),i1=1,&
!num_norbsym1(active_orbs(i)))
!enddo
!print*,'atm_nb_sym',(atm_nb_sym(i),i=1,nactorb)
!print*,'active_orb',(active_orbs(i),i=1,nactorb)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! if we take all 'pi' orbitals have non-zero overlapping!!!!

do i=1,nactorb
if(atm_nb_sym(i).lt.4)atm_nb_sym(i)=1
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,atom
print*,'atoset_geo',(atoset(i,j),j=1,atn(i))
enddo
!endif
!print*,'active_orbs',nactorb,(active_orbs(i),i=1,nactorb)
!print*,'atn',atom,(atn(i),i=1,atom)
print*,'active_atom',atom,(active_atoms(i),i=1,atom)
!print*,'atm_nb_orbs**',(atm_nb_orbs(i),i=1,nactorb)
990 format(a,x,I2,x,a,x,I2,x,a)

print*,'tot_atom',tot_atom,atom
j=0
do ij=1,tot_atom-1
!print*,'i',ij
do ii1=ij+1,tot_atom
!print*,'i1',ii1
ddist=sqrt((coordx(ij)-coordx(ii1))**2.0+(coordy(ij)-coordy(ii1))**2.0+&
(coordz(ij)-coordz(ii1))**2.0)
dist_mat(ij,ii1)=ddist
bond_dist=(at_covrad(all_at_num(ij))+at_covrad(all_at_num(ii1)))/100.0
!print*,'ddist,bond_dist',ddist,bond_dist
if(ddist.le.bond_dist)then
k1=0
do i2=1,atom
if(ij.eq.active_atoms(i2))k1=k1+1
if(ii1.eq.active_atoms(i2))k1=k1+1
enddo
if(k1.ne.2)then
dist_act_rel_mat(ij,ii1)=0.0
else
dist_act_rel_mat(ij,ii1)=1.0
endif
dist_rel_mat(ij,ii1)=1.0
else
k1=0
do i2=1,atom
if(ij.eq.active_atoms(i2))k1=k1+1
if(ii1.eq.active_atoms(i2))k1=k1+1
enddo
if(k1.ne.2)then
dist_act_rel_mat(ij,ii1)=0.0
else
dist_act_rel_mat(ij,ii1)=ddist/bond_dist
endif
dist_rel_mat(ij,ii1)=ddist/bond_dist
endif
enddo
enddo


do i=1,tot_atom
dist_act_rel_mat(i,i)=0.0
dist_rel_mat(i,i)=0.0
dist_mat(i,i)=0.0
enddo

do i=1,tot_atom-1
do i1=i+1,tot_atom
dist_act_rel_mat(i1,i)=dist_act_rel_mat(i,i1)
dist_rel_mat(i1,i)=dist_rel_mat(i,i1)
dist_mat(i1,i)=dist_mat(i,i1)
enddo
enddo

!do i=1,tot_atom
!print*,'dist_mat',(dist_mat(i,i1),i1=1,tot_atom)
!enddo
!
!do i=1,tot_atom
!print*,'dist_rel_mat',(dist_rel_mat(i,i1),i1=1,tot_atom)
!enddo
!
!do i=1,tot_atom
!print*,'dist_act_rel_mat',(dist_act_rel_mat(i,i1),i1=1,tot_atom)
!enddo
!do i=1,k
!print*,'read_info:nnmat_act',(nnmat_act(i,i1),i1=1,2)
!enddo
!do i=1,k2
!print*,'read_info:nnmat_inact',(nnmat_inact(i,i1),i1=1,2)
!enddo

!!!!!!!!!!!! below active_atoms() moved to the regular order as
!given in INFO ... for dist matrix formation it was (may be) necessary
! to change it in the order of active orbitals.

j=0

413 k=10000
do i=1,atom
do i1=1,j
if(active_atoms_sort(i1).eq.active_atoms(i))goto 412
enddo
if(k.gt.active_atoms(i))k=active_atoms(i)
412 enddo
j=j+1
active_atoms_sort(j)=k
print*,'k',k
if(j.ne.atom)goto 413

do i=1,atom
active_atoms(i)=active_atoms_sort(i)
enddo

print*,'active_atoms_sort',(active_atoms(i),i=1,atom)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

k=0
do i=1,tot_atom
least=1000.0
do i1=1,tot_atom
if(dist_act_rel_mat(i,i1).lt.least.and.dist_act_rel_mat(i,i1).ne.0.0)least=dist_act_rel_mat(i,i1)
enddo
do i1=1,tot_atom
k1=0
do i2=1,atom
if(i.eq.active_atoms(i2))k1=k1+1
if(i1.eq.active_atoms(i2))k1=k1+1
enddo
if(k1.eq.2.and.dist_act_rel_mat(i,i1).eq.least)then
k=k+1
nnmat_act(k,1)=i
nnmat_act(k,2)=i1
endif
enddo
enddo

do i=1,k
nnact(i)=nnmat_act(i,1)**2.0+nnmat_act(i,2)**2.0
enddo

do i=1,k
print*,'nnact',nnact(i),nnmat_act(i,1),nnmat_act(i,2)
enddo

j1=0
do i=1,k
do i1=1,j1
if(nnact(i).eq.nnact(nnact_1(i1)))goto 233
enddo
j1=j1+1
nnact_1(j1)=i
233 enddo
nnnatom=j1
do i=1,nnnatom
do j=1,2
nnat_bond(i,j)=nnmat_act(nnact_1(i),j)
!print*,'nnat_bond',nnat_bond(i,j)
enddo
enddo


k2=0
do i=1,tot_atom
least=1000.0
do i1=1,tot_atom
if(dist_rel_mat(i,i1).lt.least.and.dist_rel_mat(i,i1).ne.0.0)least=dist_rel_mat(i,i1)
enddo
do i1=1,tot_atom
k1=0
do i2=1,atom
if(i.eq.active_atoms(i2))k1=k1+1
if(i1.eq.active_atoms(i2))k1=k1+1
enddo

if(k1.ne.2.and.dist_rel_mat(i,i1).eq.least)then
k2=k2+1
nnmat_inact(k2,1)=i
nnmat_inact(k2,2)=i1
endif
enddo
enddo
do i=1,100
nnact(i)=0
nnact_1(i)=0
enddo

do i=1,k2
nnact(i)=nnmat_inact(i,1)**2.0+nnmat_inact(i,2)**2.0
enddo

j1=0
do i=1,k2
do i1=1,j1
if(nnact(i).eq.nnact(nnact_1(i1)))goto 234
enddo
j1=j1+1
nnact_1(j1)=i
234 enddo
nnatominact=j1
do i=1,nnatominact
do j=1,2
nnat_bond_inact(i,j)=nnmat_inact(nnact_1(i),j)
enddo
enddo

!!do i=1,tot_atom
!!print*,'dist_mat',(dist_mat(i,i1),i1=1,tot_atom)
!!enddo
!!
!!do i=1,tot_atom
!!print*,'dist_rel_mat',(dist_rel_mat(i,i1),i1=1,tot_atom)
!!enddo
!!
!!do i=1,tot_atom
!!print*,'dist_act_rel_mat',(dist_act_rel_mat(i,i1),i1=1,tot_atom)
!!enddo
!!do i=1,k
!!print*,'read_info:nnmat_act',(nnmat_act(i,i1),i1=1,2)
!!enddo
!!do i=1,k2
!!print*,'read_info:nnmat_inact',(nnmat_inact(i,i1),i1=1,2)
!!enddo


do i=1,nnnatom
print*,'read_info:nnat_bond',nnnatom,(nnat_bond(i,i1),i1=1,2)
enddo
print*,'nnatominact',nnatominact
do i=1,nnatominact
print*,'read_info:nnat_bond_inact',nnatominact,(nnat_bond_inact(i,i1),i1=1,2)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!print*,'tot_atom:geo',tot_atom
n=0
n1=0
do i=1,tot_atom
if(atoset(i,1).eq.0)goto 189
!print*,'**********',atoset(i,1),i,symat(i)
n1=n1+1
do i3=1,88
if(symat(i).eq.at_list(i3))then
!print*,'symat(i).eq.at_list(i3)',symat(i),at_list(i3)
n=n+1
actv_atom(n)=i
at_num(n)=i3
if(i3.le.54)val_st(n)=valence_state(i3)
endif
enddo
189 enddo
if(n1.ne.n)then
n=0
do i=1,tot_atom
if(atoset(i,1).eq.0)goto 190
!print*,'**********',atoset(i,1),i,symat(i)
do i3=1,88
if(symat(i).eq.at_list_bold(i3))then
!print*,'symat(i).eq.at_list(i3)',symat(i),at_list_bold(i3)
n=n+1
actv_atom(n)=i
at_num(n)=i3
if(i3.le.54)val_st(n)=valence_state(i3)
endif
enddo
190 enddo

endif

do i=1,n
print*,'at_num:geo',at_num(i),actv_atom(i),val_st(i)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Routine to produce the the inactive bias to active atooms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i1=1,20
biasval(i1)=0.0
enddo

do i=1,atom
  do i1=1,100
    biasmat(i1)=0
  enddo
k=0
k=k+1
biasmat(k)=active_atoms(i)
  do i1=1,nnatominact
    do j=1,2
      if(active_atoms(i).eq.nnat_bond_inact(i1,j))goto 235
    enddo
  goto 239
235 if(j.eq.1)then
     k=k+1
     biasmat(k)=nnat_bond_inact(i1,2)
     biasval(i)=biasval(i)+all_at_num(nnat_bond_inact(i1,2))+dist_mat(active_atoms(i),nnat_bond_inact(i1,2))
    else
     k=k+1
     biasmat(k)=nnat_bond_inact(i1,1)
     biasval(i)=biasval(i)+all_at_num(nnat_bond_inact(i1,1))+dist_mat(active_atoms(i),nnat_bond_inact(i1,1))
    endif
239 enddo
ncatm=k
nncatm=0
237 do i2=1+nncatm,ncatm
      do i1=1,nnatominact
        do j=1,2
          if(biasmat(i2).eq.nnat_bond_inact(i1,j))goto 236
        enddo
      goto 238
236 if(j.eq.1)then
        do i3=1,k
          if(biasmat(i3).eq.nnat_bond_inact(i1,2))goto 238
        enddo
        k=k+1
        biasmat(k)=nnat_bond_inact(i1,2)
     biasval(i)=biasval(i)+all_at_num(nnat_bond_inact(i1,2))+dist_mat(active_atoms(i),nnat_bond_inact(i1,2))
     else
         do i3=1,k
           if(biasmat(i3).eq.nnat_bond_inact(i1,1))goto 238
         enddo
         k=k+1
         biasmat(k)=nnat_bond_inact(i1,1)
     biasval(i)=biasval(i)+all_at_num(nnat_bond_inact(i1,1))+dist_mat(active_atoms(i),nnat_bond_inact(i1,1))
     endif
238 enddo
  enddo
nncatm=ncatm
ncatm=k
if(nncatm.ne.ncatm)goto 237
enddo

do i=1,atom
print*,'biasval_active',biasval(i)
enddo
do i=1,nnnatom
print*,'nnat_bond_ild',(nnat_bond(i,j),j=1,2)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Routine to take care of the Islands of active atooms !!!!
!!!!!!!!! Islands means the cluster of atoms !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i1=1,50
do i=1,20
islands(i,i1)=0
enddo
enddo

i4=0
do i=1,atom
if(i.eq.1)goto 246
do i1=1,i4
do j=1,nisland(i1)
if(active_atoms(i).eq.islands(i1,j))goto 245
enddo
enddo
246 k=0
k=k+1
i4=i4+1
islands(i4,k)=active_atoms(i)
do i1=1,nnnatom
do j=1,2
if(active_atoms(i).eq.nnat_bond(i1,j))goto 240
enddo
goto 241
240 if(j.eq.1)then
k=k+1
islands(i4,k)=nnat_bond(i1,2)
else
k=k+1
islands(i4,k)=nnat_bond(i1,1)
endif
241 enddo
ncatm=k
nncatm=0
242 do i2=1+nncatm,ncatm
do i1=1,nnnatom
do j=1,2
if(islands(i4,i2).eq.nnat_bond(i1,j))goto 243
enddo
goto 244
243 if(j.eq.1)then
do i3=1,k
if(islands(i4,i3).eq.nnat_bond(i1,2))goto 244
enddo
k=k+1
islands(i4,k)=nnat_bond(i1,2)
else
do i3=1,k
if(islands(i4,i3).eq.nnat_bond(i1,1))goto 244
enddo
k=k+1
islands(i4,k)=nnat_bond(i1,1)
endif
244 enddo
enddo
nncatm=ncatm
ncatm=k
!print*,'ncatm,nncatm',ncatm,nncatm
if(nncatm.ne.ncatm)goto 242
nisland(i4)=k
!print*,'nisland(i4)',nisland(i4)
245 enddo
nland=i4
!print*,'nland',nland
!do i=1,nland
!print*,'islands',(islands(i,i1),i1=1,nisland(i))
!enddo




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!lowest distance between two islands
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(nland.gt.1)then
j1=nnnatom
do i=1,nland
island_mat(i,i)=0.0
enddo
do i=1,nland-1
do j=i+1,nland
least=100.0
do i1=1,nisland(i)
do i2=1,nisland(j)
if(dist_act_rel_mat(islands(i,i1),islands(j,i2)).lt.least.and.dist_act_rel_mat&
(islands(i,i1),islands(j,i2)).ne.0.0)least=dist_act_rel_mat(islands(i,i1),islands(j,i2))
enddo
enddo
!print*,'island:least',least
island_mat(i,j)=least
island_mat(j,i)=least
enddo
enddo

do i=1,nland
print*,'island_mat(i,j)',(island_mat(i,j),j=1,nland)
enddo

do i=1,nland
least=100.0
do j=1,nland
if(island_mat(i,j).lt.least.and.island_mat(i,j).ne.0.0)then
least=island_mat(i,j)
k1=i
k2=j
island_mat(j,i)=0.0
endif
enddo

dO i1=1,nisland(k1)
do i2=1,nisland(k2)
if(dist_act_rel_mat(islands(k1,i1),islands(k2,i2)).eq.least)then
j1=j1+1
nnat_bond(j1,1)=islands(k1,i1)
nnat_bond(j1,2)=islands(k2,i2)

endif
enddo
enddo
enddo
nnnatom=j1
endif

do i=1,nnnatom
print*,'nnat_bond',(nnat_bond(i,j),j=1,2)
enddo

do i=1,atom
kk=0.0
do j=1,nnnatom
do k=1,2
if(active_atoms(i).eq.nnat_bond(j,k))then
if(k.eq.1)l=2
if(k.eq.2)l=1
kk=kk+dist_mat(active_atoms(i),nnat_bond(j,l))
endif
enddo
enddo
dist_nnat(i)=kk
print*,'dist_nnat',i,kk
enddo

do i=1,nactorb
do j=1,nactorb
if(atm_nb_orbs(i).eq.atm_nb_orbs(j))then
if(atm_nb_sym(i).eq.4.and.atm_nb_sym(j).eq.4)then
sig_sym_flg=1
goto 550
endif
endif
enddo
enddo

print*,'exit geocal'
550 return
end subroutine geocal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Keycheck(line)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
integer::i,j,k,l,nkey
character(len=200)::line,line1
character(len=10)::Keywd(20),lower(26)
data Keywd/'nset','mout','loose','tight','check','sym','chinst','ovlp','btos',&
'stob','qual','iab','nnb','sbb','udr','udb','prad','lpst','imbd','imps'/
data lower/'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q',&
'r','s','t','u','v','w','x','y','z'/

!print*,line
line1=''
nkey=20
l=len(trim(line))
k=0
do j=1,l+1
!print*,'line(j:j)',j,'>',line(j:j)
do i=1,26
if(line(j:j).eq.lower(i))then
k=k+1
line1(k:k)=lower(i)
!print*,'line1',k,line1(k:k)
goto 100
endif
enddo
!print*,'line1',line1,line,len(trim(line1))
if(len(trim(line1)).gt.0)then
do i=1,nkey
if(trim(line1).eq.Keywd(i))then
line1=''
k=0
goto 100
endif
enddo
write(*,111),'UNKNOWN KEYWORD:',trim(line1)
write(*,112),'Please choose the keywords from ',(keywd(i),i=1,20)
stop
else
goto 100
endif
100 enddo

111 format(a,2x,a)
112 format(a,2x,20a)


return
end subroutine Keycheck
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function lowercase(line)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
integer::i,j,k,l
character(len=200)::line,lowercase
character(len=2)::line1,lower(26),upper(26)
data upper/'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q'&
,'R','S','T','U','V','W','X','Y','Z'/
data lower/'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q',&
'r','s','t','u','v','w','x','y','z'/

!print*,'line1',line
l=len(line)

do j=1,l
line1=line(j:j)
do i=1,26
if(line1.eq.upper(i))then
line(j:j)=lower(i)
goto 100
endif
enddo
100 enddo

!print*,'line2',line
lowercase=line
return
end function lowercase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use commondat
use commondat1
implicit none

integer::i,MDP,io,a,b,i1,j,j1,i2,i3,d,i5,i6,i4,i7,st_num1(100),st_num2(100),num(500),m5,m6,k,kk,p,q,&
l,ll,stn,k1,k2,k3,k4,qult1(100),nbond,nel,n,nset(10),pent_set(5),info_val(4),argnum,num_2,ovlp,runn
character(len=4)::abc,at_orb(100),im_bnd(100),ato,str,at_bas_sym(100),strc,rumr,info_set(4)
character(len=200)::charst(1000)
character(len=2)::prio
character(len=200)::line34,line35,line5,lowercase
character(len=5)::pent,stn1,stn2,stn3
!character(len=15)::sttr,charst(100)
character(len=200)::sttr,Seq(20),aa
real*8::num_1(500)
logical :: fileexists
character(len=35)::inputfilename
integer::atsymset(20,20),nsym,syn(50),at_sym(50),MDP1

common/ats/atsymset,nsym,syn,at_sym
print*,'enter input_1'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! default values !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

symtype='loose'
set_order=0
iab_length=0.0
noq0=1000
noq1=1000
noq2=1000
noq3=1000
symm=0
ovval=1.000
ovopt=0
nlpset=0
nfset=0
runn=0
mset=2
qflg=0
!mult=1
ovlp_int=0
  if(input_flg.eq.1)flg1=1
 if(input_flg.eq.1)nsym=1
if(input_flg.eq.1)nstrt=0
 if(input_flg.eq.1)flgst=1
 if(input_flg.eq.1)niach=0
 if(input_flg.eq.1)niabd=0
 if(input_flg.eq.1)nialp=0
!do i=1,10
!nset(i)=100
!enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MDP=0
i6=0
i7=0
n=0

if(input_flg.eq.1)then
open(unit=21,file='input.dat',status='unknown')


do
read(21,'(a)',iostat=io)

if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(21)

do i=1,MDP
read(21,'(a)')charst(i)
print*,charst(i)
enddo

rewind(21)
endif
!! do
!!print*,'sourav'
!!      read( 21, *, iostat = IO )sttr    ! reads to EOF -- GOOD!!
!!      if ( IO < 0 ) exit               ! if EOF is reached, exit do
!!      i = i + 1
!!!      NumRec = i                       ! used later for total no. of records
!!!      allocate( Seq(I) )
!!      Seq(I) = sttr
!!      print*, I, Seq(I)
!!      X = Len_Trim( Seq(I) )           ! length of individual sequence
!!      print*, 'Sequence size: ', X
!!      print*
!!    end do
!!stop


if(input_flg.eq.0)then

call getarg(1,inputfilename)
INQUIRE(FILE=TRIM(inputfilename),EXIST=fileexists)
IF (fileexists) THEN
open(unit=21,file=TRIM(inputfilename),status='old')
ELSE
write(*,*)'SORRY This input file does not exist or you may not provide the &
filename at all'
stop
ENDIF

do
read(21,'(a)',iostat=io)

if(io.ne.0)exit
MDP=MDP+1
enddo
rewind(21)

do i=1,MDP
read(21,'(a)')charst(i)
if(index(lowercase(charst(i)),'$end').ne.0)then
MDP1=i
endif
enddo

rewind(21)
endif




do i=MDP1+1,MDP

kk=0
q=0
do j=1,100
st_num1(j)=0
st_num2(j)=0
enddo
do j=1,500
num(j)=0
enddo

!ll=len(trim(charst(i)))
!do k=1,ll+1
!line5=charst(i)(k:k)
!if(line5.eq.'')then
!sttr=trim(charst(i)(1:k-1))
!goto 111
!endif
!enddo

if(index(charst(i),'#').eq.1)goto 789
if(index(charst(i),'#').gt.1)then
ll=index(charst(i),'#')
do k=1,ll
line5=charst(i)(k:k)
if(line5.eq.'')then
!print*,'charst(i)(1:k-1)',charst(i)(1:k-1)
if(trim(charst(i)(1:k-1)).ne.'')then
charst(i)=charst(i)(1:ll-1)
goto 111
endif
endif
enddo
endif


!111 if (sttr.eq.'info'.and.input_flg.eq.1)then
111 charst(i)=lowercase(charst(i))
call Keycheck(charst(i))

if (index(charst(i),'info').ne.0.and.input_flg.eq.1)then
do j=1,i-1
read(21,*)
enddo
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
q=0
kk=0
p=0
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 133
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 143
133 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
143 enddo
stn=0
do k=1,500
num(k)=0
enddo
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
!if(fg10.eq.1)goto 189
if(line34(a:b).ne.'info')then
if(line34(a:b).eq.'nao')then 
stn=stn+1
info_set(stn)='nao'
goto 169
endif
if(line34(a:b).eq.'nae')then 
stn=stn+1
info_set(stn)='nae'
goto 169
endif
if(line34(a:b).eq.'mult')then 
stn=stn+1
info_set(stn)='mult'
goto 169
endif
if(line34(a:b).eq.'niao')then 
stn=stn+1
info_set(stn)='niao'
goto 169
endif
if(stn.ne.0)then
read(line34(a:b),'(I10)')num(stn)
!print*,num(stn),stn
info_val(stn)=num(stn)
!print*,'info_val(stn)',info_val(stn),stn,info_set(stn)
endif
endif
endif
169 enddo
!if(stn.ne.4)then
do k=1,stn
if(info_set(k).eq.'nao')then
nao=info_val(k)
goto 171
endif
enddo
write(*,*)'number of active orbitals "nao" is not given'
stop 
171 do k=1,stn
if(info_set(k).eq.'nae')then
nae=info_val(k)
goto 172
endif
enddo
write(*,*)'number of active orbitals "nae" is not given'
stop 
172 do k=1,stn
if(info_set(k).eq.'mult')then
mult=info_val(k)
goto 173
endif
enddo
write(*,*)'number of active orbitals "mult" is not given'
stop 
173 do k=1,stn
if(info_set(k).eq.'niao')then
niao=info_val(k)
goto 174
endif
enddo
write(*,*)'number of active orbitals "niao" is not given'
stop 
174 rewind(21)
endif 



!if (sttr.eq.'ovlp')then
if (index(charst(i),'ovlp').ne.0)then
do j=1,i-1
read(21,*)
enddo
!read(21,*)abc,symm

read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
p=0
q=0
kk=0
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 53
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 54
53 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
54 enddo
stn=0
do k=1,500
num_1(k)=0
enddo
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
!if(fg10.eq.1)goto 389
if(line34(a:b).eq.'ovlp')ovlp=1
if(line34(a:b).ne.'ovlp')then
if(line34(a:b).ne.'')then
stn=stn+1
!if(line34(a:b).ne.'-')then
if(stn.eq.2)then
read(line34(a:b),*)num_1(stn)
ovval=num_1(stn)
endif
if(stn.eq.1)then
read(line34(a:b),'(I10)')num_2
ovlp=num_2
endif
!print*,num(stn),stn
!endif
endif
endif
endif
enddo
print*,'ovval',ovval,ovlp

rewind(21)
endif

!if (sttr.eq.'sym')then
if (index(charst(i),'sym').ne.0)then
symm=1
!do j=1,i-1
!read(21,*)
!enddo

!read(21,'(a)')line34
if(index(charst(i),'loose').ne.0)then
symtype='loose'
endif
if(index(charst(i),'tight').ne.0)then
symtype='tight'
endif
if(index(charst(i),'check').ne.0)then
symtype='check'
endif
!stn1='qual'
!stn2='s2b'
!stn3='b2s'
if(index(charst(i),'qual').ne.0)then
set_order=0
endif
if(index(charst(i),'stob').ne.0)then
set_order=1
endif
if(index(charst(i),'btos').ne.0)then
set_order=2
endif

!print*,'set_order',set_order
rewind(21)
endif

if(index(charst(i),'inactive_bonds').ne.0.and.symm.eq.1.and.input_flg.eq.1)then
do j=1,i-1
read(21,*)
enddo
read(21,'(a)')line34
print*,line34
line35=line34
ll=len(trim(line34))
p=0
q=0
kk=0
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 2739
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 2740
2739 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
2740 enddo
stn=0
do k=1,500
num(k)=0
enddo
m5=0
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'inactive_bonds')then
stn=stn+1
if(line34(a:b).ne.'-')then
read(line34(a:b),'(I10)')num(stn)
m5=m5+1
iabd(m5)=num(stn)
endif
endif
endif
enddo
niabd=m5
!print*,'iabd(m5)',(iabd(k),k=1,m5),niabd
rewind(21)
endif


!if(sttr.eq.'inactive_lnps'.and.symm.eq.1.and.input_flg.eq.1)then
if(index(charst(i),'inactive_lnps').ne.0.and.symm.eq.1.and.input_flg.eq.1)then

do j=1,i-1
read(21,*)
enddo
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
p=0
q=0
kk=0
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 3739
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 3740
3739 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
3740 enddo
stn=0
do k=1,500
num(k)=0
enddo
m5=0
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'inactive_lnps')then
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
m5=m5+1
ialp(m5)=num(stn)
endif
endif
enddo
nialp=m5
rewind(21)
endif

!if(sttr.eq.'inactive_chrg'.and.symm.eq.1.and.input_flg.eq.1)then
if(index(charst(i),'inactive_chrg').ne.0.and.symm.eq.1.and.input_flg.eq.1)then

do j=1,i-1
read(21,*)
enddo
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
p=0
q=0
kk=0
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 339
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 340
339 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
340 enddo
stn=0
do k=1,500
num(k)=0
enddo
m5=0
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'inactive_chrg')then
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
m5=m5+1
iach(m5)=num(stn)
endif
endif
enddo
niach=m5
rewind(21)
endif


!if (sttr.eq.'frag_cntr')then
if (index(charst(i),'frag_cntr').ne.0)then
do j=1,i-1
read(21,*)
enddo
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 3841
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 3842
3841 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
3842 enddo
stn=0
do k=1,500
num(k)=0
enddo
m5=0
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'frag_cntr')then
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
m5=m5+1
frag_cntr(m5)=num(stn)
endif
endif
enddo
num_frag_cntr=m5
!print*,'frag_cntrs =',(frag_cntr(k),k=1,num_frag_cntr),'num_frag_cntr =',num_frag_cntr
rewind(21)
endif 



if(index(charst(i),'mout').ne.0)then
do j=1,i-1
read(21,*)
enddo
Read(21,*),aa,mset
if(mset.eq.0)mset=1
mset=mset+1
rewind(21)
endif



!if (sttr.eq.'nset')then
if(index(charst(i),'nset').ne.0)then

do j=1,i-1
read(21,*)
enddo
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 3839
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 3840
3839 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
3840 enddo
stn=0
do k=1,500
num(k)=0
enddo
m5=0
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'nset')then
if(line34(a:b).ne.'')then
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
m5=m5+1
nset(m5)=num(stn)
if(m5.eq.1)nfset=nset(1)
!if(nset(1).ne.2)then
if(nset(1).ne.0)then
!print*,'m555555555',m5,nset(m5)
 if(m5.eq.2)noq0=nset(2)
 if(m5.eq.3)noq1=nset(3)
 if(m5.eq.4)noq2=nset(4)
 if(m5.eq.5)noq3=nset(5)
endif
endif
endif
endif
enddo
print*,'noq0,noq1,noq2,noq3',noq0,noq1,noq2,noq3,nfset
!print*,'nset',(nset(k),k=1,5)
rewind(21)
endif 
!print*,'nfset',abc,nfset


!if (sttr.eq.'atoms'.and.input_flg.eq.1)then
if (index(charst(i),'atoms').ne.0.and.input_flg.eq.1)then
do j=1,i-1
read(21,*)
enddo
!read(21,*)abc,atom
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
kk=0
p=0
q=0
do k=1,10
pent_set(k)=0
enddo
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 263
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 265
263 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
265 enddo
stn=0
do k=1,500
num(k)=0
enddo
m5=0
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'atoms')then
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
m5=m5+1
pent_set(m5)=num(stn)
endif
endif
enddo
atom=pent_set(1)
if(atom.eq.0)then
write(*,*)'plese write the number of active atoms (atoms corresponding to &
the active orbitals) after the keyword "atoms" and give the orbital &
numbers corresponding to the each atom in the next line' 
stop
endif
!print*,'atom',atom

do i2=1,atom
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
p=0
q=0
kk=0
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 139
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 140
139 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
140 enddo
stn=0
do k=1,500
num(k)=0
enddo
n=0
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
do i3=1,88
if(line34(a:b).eq.at_list(i3))then
n=n+1
at_num(i2)=i3
goto 189
elseif(line34(a:b).eq.at_list_bold(i3))then
n=n+1
at_num(i2)=i3
!print*,'at_num**',at_num(i2)
goto 189
endif
enddo
if(n.eq.0)then
write(*,991),'Your atom',line34(a:b),'is not included in our list (atomic number 1 to &
88 are included) or you have just made a mistake'
stop
endif
991 format(a,x,a,x,a)
if(line34(a:b).ne.'-')then
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
!print*,num(stn),stn
atoset(i2,stn)=num(stn)
active_atoms(i2)=i2
endif
endif
189 enddo
atn(i2)=stn
atn(i6)=m5
!print*,'atoset_sr',i2,atn(i2),(atoset(i2,i7),i7=1,atn(i2))
!print*,'at_num',at_num(n)
enddo
rewind(21)
!do i3=1,atom
!print*,'at_num:input_1',at_num(i3)
!enddo
endif

!if (sttr.eq.'orbsym'.and.input_flg.eq.1)then
if (index(charst(i),'orbsym').ne.0.and.input_flg.eq.1)then
do j=1,i-1
read(21,*)
enddo
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
kk=0
p=0
q=0
do k=1,10
pent_set(k)=0
enddo
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 25
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 26
25 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
26 enddo
stn=0
do k=1,500
num(k)=0
enddo
m5=0
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'orbsym')then
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
m5=m5+1
pent_set(m5)=num(stn)
endif
endif
enddo
nsym=pent_set(1)
if(nsym.ge.1)then
!do i2=1,nsym
!i3=i2
!write(str,'(I1)') i3
!at_sym(i2)=trim('sym')//trim(str)
!enddo
do i4=1,nsym
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
p=0
q=0
kk=0
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 27
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 28
27 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
28 enddo
stn=0
do k=1,500
num(k)=0
enddo
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
!if(fg10.eq.1)goto 289
!print*,'at_sym(i4)',i4,at_sym(i4)
if(line34(a:b).eq.'pi')at_sym(i4)=1
if(line34(a:b).eq.'sig')at_sym(i4)=2
if(line34(a:b).ne.'pi')then
if(line34(a:b).ne.'sig')then
stn=stn+1
if(line34(a:b).ne.'-')then
!print*,'line34(a:b)',line34(a:b)
read(line34(a:b),'(I10)')num(stn)
!print*,num(stn),stn
endif
endif
endif
endif
enddo
m5=0
do k=1,stn
if(num(k).eq.0)then
do m6=num(k-1)+1,num(k+1)-1
m5=m5+1
atsymset(i4,m5)=m6
enddo
endif
if(num(k).ne.0)then
m5=m5+1
atsymset(i4,m5)=num(k)
endif
enddo
syn(i4)=m5
!print*,noqult,(qult(k),k=1,noqult)
!read(21,*)abc,a,(atsymset(i7,i5),i5=1,syn(i7))
!print*,'syn,symset',syn(i4),at_sym(i4),(atsymset(i4,i5),i5=1,syn(i4))
enddo
endif
rewind(21)
endif


!if (sttr.eq.'at_basis_sym'.and.input_flg.eq.1)then
if (index(charst(i),'at_basis_sym').ne.0.and.input_flg.eq.1)then
do j=1,i-1
read(21,*)
enddo
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
kk=0
p=0
q=0
do k=1,10
pent_set(k)=0
enddo
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 253
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 255
253 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
255 enddo
stn=0
do k=1,500
num(k)=0
enddo
m5=0
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'at_basis_sym')then
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
m5=m5+1
pent_set(m5)=num(stn)
endif
endif
enddo
nabsym=pent_set(1)

!print*,'nabset',nabsym
if(nabsym.ge.1)then
do i2=1,nabsym
i3=i2
write(str,'(I1)') i3
at_bas_sym(i2)=trim('sym')//trim(str)
enddo
do i4=1,nabsym
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
p=0
q=0
kk=0
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 239
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 240
239 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
240 enddo
stn=0
do k=1,500
num(k)=0
enddo
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
!if(fg10.eq.1)goto 289
!print*,'at_bas_sym(i4)',i4,at_bas_sym(i4)
if(line34(a:b).ne.at_bas_sym(i4))then
stn=stn+1
if(line34(a:b).ne.'-')then
!print*,'line34(a:b)',line34(a:b)
read(line34(a:b),'(I10)')num(stn)
!print*,num(stn),stn
endif
endif
endif
enddo
m5=0
do k=1,stn
if(num(k).eq.0)then
do m6=num(k-1)+1,num(k+1)-1
m5=m5+1
at_ab_symset(i4,m5)=m6
enddo
endif
if(num(k).ne.0)then
m5=m5+1
at_ab_symset(i4,m5)=num(k)
endif
enddo
absyn(i4)=m5
!print*,noqult,(qult(k),k=1,noqult)
!read(21,*)abc,a,(atsymset(i7,i5),i5=1,syn(i7))
!print*,'syn,at_ab_symset',absyn(i4),(at_ab_symset(i4,i5),i5=1,absyn(i4))
enddo
endif
rewind(21)
endif


!if (sttr.eq.'strc'.and.input_flg.eq.1)then
if (index(charst(i),'strc').ne.0.and.input_flg.eq.1)then
!print*,'sttr',sttr,i
do j=1,i-1
read(21,*)
enddo
read(21,'(a)')line34
line35=line34
!print*,line35
ll=len(trim(line34))
kk=0
p=0
q=0
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 539
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 540
539 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
540 enddo
stn=0
do k=1,500
num(k)=0
enddo
m5=0
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'strc')then
abc=line34(a:b)
endif
endif
enddo


rewind(21)
if(abc.eq.'full') flgst=1
if(abc.eq.'all') flgst=1
if(abc.eq.'covt') flgst=2
if(abc.eq.'cov') flgst=2
if(abc.eq.'ionc') flgst=3
if(abc.eq.'ion') flgst=3
!print*,flgst
endif 



!if (sttr.eq.'imbd')then
if (index(charst(i),'imbd').ne.0)then
do j=1,i-1
read(21,*)
enddo
!read(21,*)abc,imbd
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
kk=0
p=0
q=0
do k=1,10
pent_set(k)=0
enddo
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 4939
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 4940
4939 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
4940 enddo
stn=0
do k=1,500
num(k)=0
enddo
m5=0
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'imbd')then
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
m5=m5+1
pent_set(m5)=num(stn)
endif
endif
enddo
imbd=pent_set(1)
if(imbd.ne.0)then

read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
kk=0
q=0
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 1739
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 1740
1739 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
1740 enddo
stn=0
do k=1,500
num(k)=0
enddo
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'-')then
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
endif
endif
enddo
m5=0
do k=1,stn
if(num(k).ne.0)then
m5=m5+1
main_bond(m5)=num(k)
!print*,'mian_bond',main_bond(m5)
endif
enddo
nmbond=m5/2
!print*,'nmbond',nmbond
m5=0
i1=0

endif
rewind(21)
endif


!if(sttr.eq.'nnat'.and.input_flg.eq.1)then
if(index(charst(i),'nnat').ne.0.and.input_flg.eq.1)then
do j=1,i-1
read(21,*)
enddo
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 739
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 740
739 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
740 enddo
stn=0
do k=1,500
num(k)=0
enddo
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
!if(fg10.eq.1)goto 389
if(line34(a:b).ne.'nnat')then
stn=stn+1
if(line34(a:b).ne.'-')then
read(line34(a:b),'(I10)')num(stn)
!print*,num(stn),stn
endif
endif
endif
enddo
m5=0
do k=1,stn
if(num(k).ne.0)then
m5=m5+1
pbond(m5)=num(k)
endif
enddo
nbond=m5
m5=0
i1=0
do k=1,nbond,2
i2=0
i1=i1+1
do m5=k,k+1
i2=i2+1
nnat_bond(i1,i2)=pbond(m5)
enddo
enddo
nnnatom=nbond/2
!print*,'nbond',nbond
do i1=1,nnnatom
!print*,'nnbond',nnnatom,i1,(nnat_bond(i1,k),k=1,2)
enddo
rewind(21)
endif


!if(sttr.eq.'prad')then
if(index(charst(i),'prad').ne.0)then
do j=1,i-1
read(21,*)
enddo
!read(21,*)abc,prad
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
p=0
q=0
kk=0
do k=1,10
pent_set(k)=0
enddo
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 6939
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 6940
6939 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
6940 enddo
stn=0
do k=1,500
num(k)=0
enddo
m5=0
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'prad')then
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
m5=m5+1
pent_set(m5)=num(stn)
!print*,num(stn),stn
endif
endif
enddo
prad=pent_set(1)
if(prad.ne.0)then
j1=0
do j=i,i+prad-1
do k=1,100
st_num1(k)=0
st_num2(k)=0
enddo
kk=0
q=0
p=0
j1=j1+1
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 839
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 840
839 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
840 enddo
stn=0
do k=1,500
num(k)=0
enddo
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll+1)then
 a=st_num1(k)
 b=st_num2(k+1)
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
endif
enddo
m5=0
do k=1,stn
if(num(k).ne.0)then
m5=m5+1
prio_rad(j1,m5)=num(k)
endif
enddo
norad(j1)=m5
enddo
endif
rewind(21)
endif

!if(sttr.eq.'lpst')then
if(index(charst(i),'lpst').ne.0)then
do j=1,i-1
read(21,*)
enddo

read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
do k=1,10
pent_set(k)=0
enddo
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 5939
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 5940
5939 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
5940 enddo
stn=0
do k=1,500
num(k)=0
enddo
m5=0
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'lpst')then
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
m5=m5+1
pent_set(m5)=num(stn)
endif
endif
enddo
nlpset=pent_set(1)
if(nlpset.ne.0)then



do j=1,nlpset
kk=1
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 939
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 940
939 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
940 enddo
stn=0
do k=1,500
num(k)=0
enddo
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'lset')then
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
endif
endif
enddo
m5=0
do k=1,stn
if(num(k).eq.0)then
do m6=num(k-1)+1,num(k+1)-1
m5=m5+1
plpair(j,m5)=num(k)
enddo
endif
if(num(k).ne.0)then
m5=m5+1
plpair(j,m5)=num(k)
endif
enddo
lp(j)=m5
enddo
endif
rewind(21)
endif


!if(sttr.eq.'imps')then
if(index(charst(i),'imps').ne.0)then
do j=1,i-1
read(21,*)
enddo


read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
kk=0
p=0
q=0
do k=1,10
pent_set(k)=0
enddo
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 439
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 440
439 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
440 enddo
stn=0
do k=1,500
num(k)=0
enddo
m5=0
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).ne.'imps')then
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
m5=m5+1
pent_set(m5)=num(stn)
endif
endif
enddo
nstrt=pent_set(1)
if(nstrt.ne.0)then
l=0
do j=1,nstrt
kk=1
p=0
q=0
read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 1139
!print*,k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 1140
1139 kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
endif
1140 enddo
stn=0
do k=1,500
num(k)=0
enddo
!print*,'kk',kk
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
stn=stn+1
read(line34(a:b),'(I10)')num(stn)
endif
enddo
m5=0
if(stn.eq.0)goto 490
l=l+1
do k=1,stn
if(num(k).eq.0)then
do m6=num(k-1)+1,num(k+1)-1
m5=m5+1
strt_struc(l,m5)=num(k)
enddo
endif
if(num(k).ne.0)then
m5=m5+1
strt_struc(l,m5)=num(k)
endif
enddo
nel=m5
if(nel.ne.nae)then
write(*,*)'please check your starting structures, may be number of active electrons are wrong'
stop
endif
enddo
endif
490 rewind(21)
print*,'imps',stn,nel,nstrt
do m6=1,nstrt
print*,(strt_struc(m6,j),j=1,nel)
enddo
endif



!if (sttr.eq.'iab'.or.sttr.eq.'sbb'.or.sttr.eq.'nnb'.or.sttr.eq.'udr'.or.sttr.eq.'udb')then
if (index(charst(i),'iab').ne.0.or.index(charst(i),'sbb').ne.0&
.or.index(charst(i),'nnb').ne.0.or.index(charst(i),'udr').ne.0&
.or.index(charst(i),'udb').ne.0)then
runn=runn+1
if(runn.eq.1)then
itb=0
syb=0
nnb=0
radical=0
mnbond=0
endif

do j=1,i-1
read(21,*)
enddo

read(21,'(a)')line34
line35=line34
ll=len(trim(line34))
q=0
kk=0
p=0
do k=1,ll+1
line5=line34(k:k)
if(line5.eq.'')then
q=q+1
if(k.eq.1)p=1
if(p.eq.1)goto 3949
!print*,'k',k
if(q.eq.1)then
kk=1
st_num1(1)=1
endif
kk=kk+1
st_num1(kk)=k+1
st_num2(kk)=k
goto 3941
3949 kk=kk+1
!print*,'k**',k
st_num1(kk)=k+1
st_num2(kk)=k
endif
3941 enddo
stn=0
do k=1,500
num(k)=0
enddo
m5=0
!print*,'kk**',kk
do k=1,kk
a=0
b=0
if(st_num2(k+1)-st_num2(k).ne.1.and.st_num2(k).le.ll)then
 a=st_num1(k)
 b=st_num2(k+1)
if(line34(a:b).eq.'iab')then
pent=line34(a:b)
stn=stn+1
goto 117
endif
if(line34(a:b).eq.'sbb')then
pent=line34(a:b)
stn=stn+1
goto 117
endif
if(line34(a:b).eq.'nnb')then
pent=line34(a:b)
stn=stn+1
goto 117
endif
if(line34(a:b).eq.'udr')then
pent=line34(a:b)
stn=stn+1
goto 117
endif
if(line34(a:b).eq.'udb')then
pent=line34(a:b)
stn=stn+1
goto 117
endif
read(line34(a:b),'(I10)')num(stn)

endif



117 enddo

if(pent.eq.'iab')itb=num(stn)
if(pent.eq.'sbb')syb=num(stn)
if(pent.eq.'nnb')nnb=num(stn)
if(pent.eq.'udr')radical=num(stn)
if(pent.eq.'udb')mnbond=num(stn)
print*,'itb,syb,nnb,radical,mnbond:input',itb,syb,nnb,radical,mnbond
endif
rewind(21)

!if (sttr.eq.'iab_len')then
if (index(charst(i),'iab_len').ne.0)then

do j=1,i-1
read(21,*)
enddo

read(21,*)line34,iab_length
rewind(21)
endif

789 enddo



!do i1=1,6
!print*,nmbond(i1),i1,(main_bond(i1,k),k=1,nmbond(i1))
!enddo
!do j=1,nlpset
!print*,nlpset,lp,(plpair(j,k),k=1,lp)
!enddo
!print*,'input_1',nao,nae,mult,niao,atn(1),atn(2)

vacorb=nao-nae
nlp=nae-nao
if(vacorb.gt.1)nlp=0
if(qult(1).eq.0)then
noqult=(1+(nae-nlast-nlp*2)/2)*(1+(nae-nlast-nlp*2)/2)
do j=1,noqult
qult(j)=j
enddo
endif

do i=1,10
pent_set(i)=0
enddo
!print*,'itb,syb,nnb,radical,mnbond:11',itb,syb,nnb,radical,mnbond
if(prad.eq.0)radical=0
if(imbd.eq.0)mnbond=0
!print*,'itb,syb,nnb,radical,mnbond:22',itb,syb,nnb,radical,mnbond
pent_set(1)=itb
pent_set(2)=syb
pent_set(3)=nnb
pent_set(4)=radical
pent_set(5)=mnbond
print*,'itb,syb,nnb,radical,mnbond:22',itb,syb,nnb,radical,mnbond

l=1
do j=1,5
if(l.lt.pent_set(j))l=pent_set(j)
enddo

!print*,'ll:pent_set',l
k=0
do j=1,l
do i=1,5
if(pent_set(i).eq.j)then
k=k+1
goto 131
!pent_set(i)=k
endif
enddo
goto 130
131 do i=1,5
if(pent_set(i).eq.j)then
pent_set(i)=k
endif
enddo
130 enddo
itb=pent_set(1)
syb=pent_set(2)
nnb=pent_set(3)
radical=pent_set(4)
mnbond =pent_set(5)
       
itb=itb+1
syb=syb+1
nnb=nnb+1
radical=radical+1
mnbond=mnbond+1
print*,'itb,syb,nnb,radical,mnbond',itb,syb,nnb,radical,mnbond
if(itb.eq.1.and.syb.eq.1.and.nnb.eq.1.and.radical.eq.1.and.mnbond.eq.1.and.flg1.ne.1)then
!itb=1
!syb=2
!nnb=3
!radical=4
!mnbond=5
qflg=1
endif
nlast=(mult-1)

if(ovlp.eq.0)ovopt=0
if(ovlp.eq.2)ovopt=1
if(ovlp.eq.1.and.nfset.eq.0)ovopt=0
if(ovlp.eq.1.and.nfset.eq.1)ovopt=1
if(ovlp.eq.1.and.nfset.eq.2)ovopt=1
if(ovlp.eq.1.and.nfset.eq.3)ovopt=1
!!!!!!!!! "vpt" is the option for user specifying overlap value "ovval". It will work
!!!!!!!!!!!!!only for vpt=1. To lock it please put any other value
vpt=1

write(7,*)'******************************************************************************************'
write(7,900)'*','active orbs =',nao,'active electrons =',nae,'multiplicity =',mult,'inactive orbs =',niao
write(7,*)'******************************************************************************************'
write(7,*)'******************************************************************************************'
900 format(a,1x,a,1x,I4,3x,a,1x,I4,3x,a,1x,I4,3x,a,1x,I4)

!print*,'itb,syb,nnb,radical,mnbond',itb,syb,nnb,radical,mnbond
!print*,'exit input_1'

return
end subroutine input
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cov_struc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig

integer::ijk,ij,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,perm_nstr,nssym,totset,&
i14,i15,i16,i17,i18,i19,i20,i21,k1,k2,k3,k4,k5,STDOUT,a,c,b,d,e,f,nnn,elporb,x,y
integer::i,j,k,l,j1,m,j2,ii,jj,iii,n1,j3,j4,j5,wig2,tnqs,lonep,bonds,allp
integer::alstr,fullcovstr(15000,20),symq(15000),sigsym(15000),tnqs_sig&
,qulsym(15000),str_quality_1(15000),str_quality_2(15000),bondq(15000),tqlty,bqlty,sqlty,symqq(15000)
real*8::factorial,symsc(15000),T1,T2
integer,dimension(:,:),allocatable::strc,fstr,strct,num
integer,dimension(:),allocatable::n,nn
!!! memory allocation !!!!!!!
x=100000
y=100

allocate(strc(x,y))
allocate(strct(x,y))
allocate(n(x))
allocate(nn(x))
allocate(num(x,y))
!!! initialization of arrays !!!
strc(x,y)=0
strct(x,y)=0
n(x)=0
nn(x)=0
num(x,y)=0
!!! ionic and covalent flags
flg_ion=0
flg_cov=1
!!!! new structures generation part starts !!!!!
vacorb=nao-nae
lonep=nae-nao
if(vacorb.gt.0)lonep=0
bonds=(nae-lonep*2-nlast)/2
d=nao
e=2
f=d-e
c=factorial(d)/(factorial(e)*factorial(f))
!!!! production of the set of bonded orbitals strats !!!!!
i4=0
do i1=1,nao-1
  do i2=i1+1,nao
    i4=i4+1
    do i3=1,2
    i5=i2
      if(i3.eq.1)i5=i1
    strc(i4,i3)=i5
    n(i4)=i1
    enddo
  enddo
enddo
totset=i4
!!!! production of the set of bonded orbitals ends !!!!!
!!!! production of the covalent bonding part start !!!!!
i=0
j=0
do i1=1,totset
j=1
nn(1)=i1
  if(j.eq.bonds)then
  i=i+1
  k4=0
    do k2=1,bonds
      do k5=1,2
      k4=k4+1
      strct(i,k4)=strc(nn(k2),k5)
      enddo
    enddo
  goto 121
  endif

do i2=i1+1,totset
j=2
l=0
  do k=1,2
    if(strc(nn(1),k).eq.n(i2))then
    l=l+1
    endif
  enddo
  if(l.eq.0)goto 142
  goto 122
142 l=0 
  do k=1,2
    do k1=1,2
      if(strc(i2,k).eq.strc(nn(1),k1))then
      l=l+1
      endif
    enddo
  enddo
  if(l.eq.0)goto 162
  goto 122
162 nn(2)=i2
if(j.eq.bonds)then
i=i+1
k4=0
do k2=1,bonds
do k5=1,2
k4=k4+1
strct(i,k4)=strc(nn(k2),k5)
enddo
enddo
!print*,'strct',(strct(i,k),k=1,k4)
goto 122
endif

do i3=i2+1,totset
!print*,'sourav',i3
j=3
l=0
do k1=1,j-1
do k=1,2
if(strc(nn(k1),k).eq.n(i3))then
l=l+1
endif
enddo
enddo
!print*,'i3,l',i2,l
if(l.eq.0)goto 143
goto 123
143 l=0 
do k2=1,j-1
 do k=1,2
do k1=1,2
if(strc(i3,k).eq.strc(nn(k2),k1))then
l=l+1
endif
enddo
enddo
enddo
if(l.eq.0)goto 163
goto 123
163 nn(3)=i3
if(j.eq.bonds)then
i=i+1
k4=0
do k2=1,bonds
do k5=1,2
k4=k4+1
!print*,strc(nn(k2),k5)
strct(i,k4)=strc(nn(k2),k5)
enddo
enddo
!print*,'strct',(strct(i,k),k=1,k4)
goto 123
endif

do i4=i3+1,totset
j=4
l=0
do k1=1,j-1
do k=1,2
if(strc(nn(k1),k).eq.n(i4))then
l=l+1
endif
enddo
enddo
if(l.eq.0)goto 144
goto 124
144 l=0 
do k2=1,j-1
 do k=1,2
do k1=1,2
if(strc(i4,k).eq.strc(nn(k2),k1))then
l=l+1
endif
enddo
enddo
enddo
if(l.eq.0)goto 164
goto 124
164 nn(4)=i4
if(j.eq.bonds)then
i=i+1
k4=0
do k2=1,bonds
do k5=1,2
k4=k4+1
strct(i,k4)=strc(nn(k2),k5)
enddo
enddo
!print*,'strct',(strct(i,k),k=1,k4)
goto 124
endif


do i5=i4+1,totset
j=5
l=0
do k1=1,j-1
do k=1,2
if(strc(nn(k1),k).eq.n(i5))then
l=l+1
endif
enddo
enddo
if(l.eq.0)goto 145
goto 125
145 l=0 
do k2=1,j-1
 do k=1,2
do k1=1,2
if(strc(i5,k).eq.strc(nn(k2),k1))then
l=l+1
endif
enddo
enddo
enddo
if(l.eq.0)goto 165
goto 125
165 nn(5)=i5
if(j.eq.bonds)then
i=i+1
k4=0
do k2=1,bonds
do k5=1,2
k4=k4+1
strct(i,k4)=strc(nn(k2),k5)
enddo
enddo
!print*,'strct',(strct(i,k),k=1,k4)
goto 125
endif

do i6=i5+1,totset
j=6
l=0
do k1=1,j-1
do k=1,2
if(strc(nn(k1),k).eq.n(i6))then
l=l+1
endif
enddo
enddo
if(l.eq.0)goto 146
goto 126
146 l=0 
do k2=1,j-1
 do k=1,2
do k1=1,2
if(strc(i6,k).eq.strc(nn(k2),k1))then
l=l+1
endif
enddo
enddo
enddo
if(l.eq.0)goto 166
goto 126
166 nn(6)=i6
!print*,'nn4',nn(6)
if(j.eq.bonds)then
!print*,'sourav2',nn(1),nn(2)
i=i+1
k4=0
do k2=1,bonds
do k5=1,2
k4=k4+1
!print*,strc(nn(k2),k5)
strct(i,k4)=strc(nn(k2),k5)
enddo
enddo
!print*,'strct',(strct(i,k),k=1,k4)
goto 126
endif

do i7=i6+1,totset
!print*,'souravi4',i7,totset
j=7
l=0
do k1=1,j-1
do k=1,2
if(strc(nn(k1),k).eq.n(i7))then
l=l+1
endif
enddo
enddo
!print*,'i3,l',i2,l
if(l.eq.0)goto 147
goto 127
147 l=0 
do k2=1,j-1
 do k=1,2
do k1=1,2
if(strc(i7,k).eq.strc(nn(k2),k1))then
l=l+1
endif
enddo
enddo
enddo
if(l.eq.0)goto 167
goto 127
167 nn(7)=i7
!print*,'nn4',nn(7)
if(j.eq.bonds)then
!print*,'sourav2',nn(1),nn(2)
i=i+1
k4=0
do k2=1,bonds
do k5=1,2
k4=k4+1
!print*,strc(nn(k2),k5)
strct(i,k4)=strc(nn(k2),k5)
enddo
enddo
!print*,'strct',(strct(i,k),k=1,k4)
goto 127
endif

127 enddo
126 enddo
125 enddo
124 enddo
123 enddo
122 enddo
121 enddo
!!!! production of the covalent bonding part of the structures ends !!!!!

alstr=i

strc(x,y)=0
n(x)=0
nn(x)=0


do i1=1,alstr
do k4=1,bonds*2
strc(i1,k4)=strct(i1,k4)
enddo
enddo
!do i1=1,i
!print*,(strc(i1,k4),k4=1,bonds*2)
!enddo

!!!!! production of the lone pairs and radical part of the structures starts !!!!!!
if(nlast.ne.0.or.lonep.ne.0)then
i=0
  do ii=1,alstr
  k5=0
  i2=0
   do i1=1,nao
     do k1=1,bonds*2
       if(i1.eq.strct(ii,k1)) goto 530
     enddo
i2=i2+1
num(ii,i2)=i1

530 enddo
allp=i2
!!!!! production of the lone pair part of the structures starts !!!!!!
if(lonep.ne.0)then
j=0
m=0
  do i1=1,allp
  j=1
  n(1)=i1
    if (j.eq.lonep) then
    i=i+1
    j4=0
      do j1=lonep,1,-1
        do j2=1,2
        j4=j4+1
        strc(i,j4)=num(ii,n(j1))
        enddo
      enddo
    j2=0
      do j3=1,bonds*2
      j2=j3+j4
      strc(i,j2)=strct(ii,j3)
      enddo
    goto 399
    endif

do i2=i1+1,allp
j=2
  if(i2.eq.n(1))goto 401
  n(2)=i2
    if (j.eq.lonep) then
    i=i+1
    j4=0
      do j1=1,lonep
        do j2=1,2
        j4=j4+1
        strc(i,j4)=num(ii,n(j1))
        enddo
      enddo
    j2=0
      do j3=1,bonds*2
      j2=j3+j4
      strc(i,j2)=strct(ii,j3)
      enddo
      goto 401
    endif

do i3=i2+1,allp
j=3
do k=1,2
if(i3.eq.n(k))goto 402
enddo
n(3)=i3
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 402
endif

do i4=i3+1,allp
j=4
do k=1,3
if(i4.eq.n(k))goto 403
enddo
n(4)=i4
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 403
endif

do i5=i4+1,allp
j=5
do k=1,4
if(i5.eq.n(k))goto 404
enddo
n(5)=i5
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 404
endif

do i6=i5+1,allp
j=6
do k=1,5
if(i6.eq.n(k))goto 405
enddo
n(6)=i6
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 405
endif

do i7=i6+1,allp
j=7
do k=1,6
if(i7.eq.n(k))goto 406
enddo
n(7)=i7
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 406
endif

do i8=i7+1,allp
j=8
do k=1,7
if(i8.eq.n(k))goto 407
enddo
n(8)=i8
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 407
endif

do i9=i8+1,allp
j=9
do k=1,8
if(i9.eq.n(k))goto 408
enddo
n(9)=i9
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 408
endif

do i10=i9+1,allp
j=10
do k=1,9
if(i10.eq.n(k))goto 409
enddo
n(10)=i10
if (j.eq.lonep) then
i=i+1
j4=0
do j1=1,lonep
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 409
endif

409 enddo
408 enddo
407 enddo
406 enddo
405 enddo
404 enddo
403 enddo
402 enddo
401 enddo
399 enddo
endif
enddo
if (lonep.ne.0) alstr=i
!!!!! production of the radical part of the structures starts !!!!!!
if(nlast.ne.0)then
  do i=1,alstr
  j4=nae-nlast
    do i1=1,nae
      do i2=1,nae-nlast
        if(strc(i,i2).eq.i1)goto 555
      enddo
      j4=j4+1
      strc(i,j4)=i1
555 enddo
enddo
endif
endif
!!!!! production of the lone pairs and radical part of the structures ends !!!!!!
deallocate(num)
990 format(30I3)
656 deallocate(strct)
!!!!!!! new structures generation part ends !!!!!!!!!!!!!!!!!!!!!!!
d=nao
e=nlp
f=nao-nlp
elporb=nao-nlp
c=factorial(d)/(factorial(e)*factorial(f))
call wigner(elporb,wig2)
deallocate(n)
deallocate(nn)

write(7,*)'  '
write(7,*)'                  covalent structures   '
write(7,*)'                  -------------------   '
write(7,*)'  '

write(7,307)' You have',alstr,' covalent structures among'
write(7,307)'them',wig2*c,' covalent structures are permissible'
write(7,*)' '
write(7,*)'      [] and {} in front of the structures specifies the quality of the structures.'
write(7,*)'1st number in [] specifies intra-atomic bond quality and 2nd number specifies ' 
write(7,*)'symmetry breaking bond quality. Where {} specifies the overall quality of the structures'
write(7,*)'*************************************************************** '
307 format(2x,a,I7,a)


do i=1,alstr
do i1=1,nae
fullcovstr(i,i1)=strc(i,i1)+niao
enddo
enddo

deallocate(strc)

if(symm.eq.1) then
if(sig_sym_flg.eq.1)call symmetry_cal_sig(nlp,fullcovstr,alstr,symsc,symq,nssym)
if(sig_sym_flg.ne.1)call symmetry_cal_pi(nlp,fullcovstr,alstr,symsc,symq,nssym)

!do ij=1,alstr
!print*,'qulsym(ij),symq(ij)',qulsym(ij),symq(ij),(fullcovstr(ij,j1),j1=1,nao-nlp+nlp*2)
!print*,'qulsym(ij),symq(ij)',qulsym(ij),symq(ij),alstr,nssym,tnqs
!enddo


nnn=0
do ij=1,tnqs
do ijk=nssym,1,-1
do j2=1,alstr
if(symq(j2).eq.ijk)then
nnn=nnn+1
if(niao.eq.0)then
write(7,900),'cov structure',nnn,')',qulsym(j2),(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
endif
if(niao.ne.0)then
write(7,901),'cov structure',nnn,')',qulsym(j2),1,':',niao,(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
endif
endif
210 enddo
write(7,*)'******************************************************'
enddo
enddo
else

do j2=1,alstr
if(niao.eq.0)then
write(7,902),'cov structure',j2,')',(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
endif
if(niao.ne.0)then
write(7,903),'cov structure',j2,')',1,':',niao,(fullcovstr(j2,j1),j1=1,nao-nlp+nlp*2)
endif

enddo
write(7,*)'******************************************************'
endif
900 format(a,I5,a,x,I3,3x,30(I5))
901 format(a,I5,a,x,I3,3x,I3,x,a,I3,30(I5))
902 format(a,I5,a,x,30(I5))
903 format(a,I5,a,x,I3,x,a,I3,30(I5))


if(nfset.eq.3)write(23,*)'********** RUMER STRUCTURES ******************'
if(nfset.eq.5)write(10,*)'**********ALL POSSIBLE RUMER STRUCTURES ******************'
if(flg1.eq.1)then
write(10,*)'********** RUMER STRUCTURES ******************'
endif
if(flg1.eq.0.and.nfset.ne.5)then
write(10,*)'******* CHEM. QUAL. STRUCTURES **************'
endif

perm_nstr=wig2*c
write(10,*),perm_nstr,' covalent structures' 
write(23,*),perm_nstr,' covalent structures' 
if(nfset.ne.5)write(10,913),'IAB','NNB','SBB','fqual'
write(23,913),'IAB','NNB','SBB','fqual'
913 format(x,a,x,a,x,a,x,a)
call str_selection(fullcovstr,nlp,alstr,perm_nstr)
!print*,'total num covstr',tncs
print*,'exit cove_struc'
return
end subroutine cov_struc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ion_struc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig
common/str/str5,nstr7

   integer::i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,perm_nstr,nssym,nnn,x,y, &
i14,i15,i16,i17,i18,i19,i20,i21,k1,k2,STDOUT,rns,jjj,m1,c,c1,d,e,f,elporb,wig2,tnqs
   integer::i,j,k,l,l2,j1,m,j2,ii,jj,iii,n1,j3,j4,j5,nnlp,nnao,nos,nofs,j6,j7,j8,ij,ijk,nstruc
   integer::alionstr,symq(15000),sigsym(15000),tnqs_sig,ion_str_p(1000,15),ald_ionstr(5000,20)&
   ,qulsym(15000),fullionstr(15000,20),str_quality_1(15000),str5(2000,20),nstr7,ll,lll,&
   str_quality_2(15000),bondq(15000),tqlty,bqlty,sqlty,str_p(1000,15),istrbn,strbn
integer::llll,str2(15000,20),str3(15000,20),str4(15000,20),space(1000),covnris(1000)&
,covris(1000),fvec(15000,1000),strn(10000),l1,m119,m20,q_fac(15000),quality_fac(15000)
integer::bonds,alstr,totset,allp,k4,k5,l3,l4,a,b,totistr
   real*8::factorial,symsc(15000)
   integer,dimension(:,:),allocatable::strc,strct,num,istrc
   integer,dimension(:),allocatable::n,nn,num1

print*,'enter ion_struc',niao
x=100000
y=100

!allocate(strc(x,y))
allocate(strct(x,y))
allocate(istrc(x,y))
allocate(num(x,y))
allocate(num1(x))
allocate(n(x))
allocate(nn(x))

!strc(x,y)=0
strct(x,y)=0
istrc(x,y)=0
n(x)=0
nn(x)=0
num(x,y)=0
num1(x)=0

write(7,*)''
write(7,*)'                              ionic structures '
write(7,*)'                              ---------------- '

flg_ion=1
flg_cov=0

symm=0

m=0
j=0
i=0
nos=0
nofs=0
nnao=0
jjj=0
!print*,'nnlp',nlp,nao,nlast
if(vacorb.gt.0)nlp=nlp+1
do nnlp=nlp+1,nlp+(nao-nlp-nlast)/2
allocate(strc(x,y))
strc(x,y)=0

bonds=(nae-nnlp*2-nlast)/2
d=nao
e=2
f=d-e
!print*,'fac',d,e,f,elporb,bonds
c=factorial(d)/(factorial(e)*factorial(f))

c=bonds
d=nae-nnlp*2
b=int(d/2)
a=2**b
e=nae-nnlp*2-2*b
totistr=factorial(d)/(factorial(e)*a*factorial(c))

!print*,'totistr',totistr
!!!! production of the set of bonded orbitals strats !!!!!

if(bonds.ne.0)then
i4=0
do i1=1,nao-1
do i2=i1+1,nao
i4=i4+1
do i3=1,2
i5=i2
if(i3.eq.1)i5=i1
strc(i4,i3)=i5
n(i4)=i1
enddo
enddo
enddo
!do j2=1,i4
!print*,'bond',j2,n(j2),(strc(j2,j1),j1=1,2)
!enddo
totset=i4
!print*,'totset',i4
!!!! production of the set of bonded orbitals ends !!!!!

!!!! production of the covalent bonding part of the structures start !!!!!
i=0
j=0
do i1=1,totset
!print*,'sourav',i1,totset
j=1
nn(1)=i1
!print*,'nn1',nn(1)
if(j.eq.bonds)then
!print*,'sourav1',j
i=i+1
k4=0
do k2=1,bonds
do k5=1,2
k4=k4+1
strct(i,k4)=strc(nn(k2),k5)
enddo
enddo
goto 121
endif

do i2=i1+1,totset
j=2
l=0
do k=1,2
if(strc(nn(1),k).eq.n(i2))then
l=l+1
endif
enddo
if(l.eq.0)goto 142
goto 122
142 l=0 
 do k=1,2
do k1=1,2
if(strc(i2,k).eq.strc(nn(1),k1))then
l=l+1
endif
enddo
enddo
if(l.eq.0)goto 162
goto 122
162 nn(2)=i2
if(j.eq.bonds)then
i=i+1
k4=0
do k2=1,bonds
do k5=1,2
k4=k4+1
!print*,strc(nn(k2),k5)
strct(i,k4)=strc(nn(k2),k5)
enddo
enddo
!print*,'strct',(strct(i,k),k=1,k4)
goto 122
endif

do i3=i2+1,totset
j=3
l=0
do k1=1,j-1
do k=1,2
if(strc(nn(k1),k).eq.n(i3))then
l=l+1
endif
enddo
enddo
if(l.eq.0)goto 143
goto 123
143 l=0 
do k2=1,j-1
 do k=1,2
do k1=1,2
if(strc(i3,k).eq.strc(nn(k2),k1))then
l=l+1
endif
enddo
enddo
enddo
if(l.eq.0)goto 163
goto 123
163 nn(3)=i3
if(j.eq.bonds)then
i=i+1
k4=0
do k2=1,bonds
do k5=1,2
k4=k4+1
strct(i,k4)=strc(nn(k2),k5)
enddo
enddo
!print*,'strct',(strct(i,k),k=1,k4)
goto 123
endif

do i4=i3+1,totset
!print*,'souravi4',i4,totset
j=4
l=0
do k1=1,j-1
do k=1,2
if(strc(nn(k1),k).eq.n(i4))then
l=l+1
endif
enddo
enddo
!print*,'i3,l',i2,l
if(l.eq.0)goto 144
goto 124
144 l=0 
do k2=1,j-1
 do k=1,2
do k1=1,2
if(strc(i4,k).eq.strc(nn(k2),k1))then
l=l+1
endif
enddo
enddo
enddo
if(l.eq.0)goto 164
goto 124
164 nn(4)=i4
!print*,'nn4',nn(4)
if(j.eq.bonds)then
!print*,'sourav2',nn(1),nn(2)
i=i+1
k4=0
do k2=1,bonds
do k5=1,2
k4=k4+1
!print*,strc(nn(k2),k5)
strct(i,k4)=strc(nn(k2),k5)
enddo
enddo
!print*,'strct',(strct(i,k),k=1,k4)
goto 124
endif


do i5=i4+1,totset
!print*,'souravi4',i5,totset
j=5
l=0
do k1=1,j-1
do k=1,2
if(strc(nn(k1),k).eq.n(i5))then
l=l+1
endif
enddo
enddo
!print*,'i3,l',i2,l
if(l.eq.0)goto 145
goto 125
145 l=0 
do k2=1,j-1
 do k=1,2
do k1=1,2
if(strc(i5,k).eq.strc(nn(k2),k1))then
l=l+1
endif
enddo
enddo
enddo
if(l.eq.0)goto 165
goto 125
165 nn(5)=i5
!print*,'nn4',nn(5)
if(j.eq.bonds)then
!print*,'sourav2',nn(1),nn(2)
i=i+1
k4=0
do k2=1,bonds
do k5=1,2
k4=k4+1
!print*,strc(nn(k2),k5)
strct(i,k4)=strc(nn(k2),k5)
enddo
enddo
!print*,'strct',(strct(i,k),k=1,k4)
goto 125
endif

do i6=i5+1,totset
!print*,'souravi4',i4,totset
j=6
l=0
do k1=1,j-1
do k=1,2
if(strc(nn(k1),k).eq.n(i6))then
l=l+1
endif
enddo
enddo
!print*,'i3,l',i2,l
if(l.eq.0)goto 146
goto 126
146 l=0 
do k2=1,j-1
 do k=1,2
do k1=1,2
if(strc(i6,k).eq.strc(nn(k2),k1))then
l=l+1
endif
enddo
enddo
enddo
if(l.eq.0)goto 166
goto 126
166 nn(6)=i6
!print*,'nn4',nn(6)
if(j.eq.bonds)then
!print*,'sourav2',nn(1),nn(2)
i=i+1
k4=0
do k2=1,bonds
do k5=1,2
k4=k4+1
!print*,strc(nn(k2),k5)
strct(i,k4)=strc(nn(k2),k5)
enddo
enddo
!print*,'strct',(strct(i,k),k=1,k4)
goto 126
endif

do i7=i6+1,totset
!print*,'souravi4',i7,totset
j=7
l=0
do k1=1,j-1
do k=1,2
if(strc(nn(k1),k).eq.n(i7))then
l=l+1
endif
enddo
enddo
!print*,'i3,l',i2,l
if(l.eq.0)goto 147
goto 127
147 l=0 
do k2=1,j-1
 do k=1,2
do k1=1,2
if(strc(i7,k).eq.strc(nn(k2),k1))then
l=l+1
endif
enddo
enddo
enddo
if(l.eq.0)goto 167
goto 127
167 nn(7)=i7
!print*,'nn4',nn(7)
if(j.eq.bonds)then
!print*,'sourav2',nn(1),nn(2)
i=i+1
k4=0
do k2=1,bonds
do k5=1,2
k4=k4+1
!print*,strc(nn(k2),k5)
strct(i,k4)=strc(nn(k2),k5)
enddo
enddo
!print*,'strct',(strct(i,k),k=1,k4)
goto 127
endif

127 enddo
126 enddo
125 enddo
124 enddo
123 enddo
122 enddo
121 enddo
!!!! production of the covalent bonding part of the structures ends !!!!!

alstr=i

strc(x,y)=0
n(x)=0
nn(x)=0


do i1=1,alstr
do k4=1,bonds*2
strc(i1,k4)=strct(i1,k4)
enddo
enddo
!do i1=1,alstr
!print*,(strc(i1,k4),k4=1,bonds*2)
!enddo
endif
if(bonds.eq.0)then
alstr=1
strc(x,y)=0
strct(x,y)=0
endif
!do i1=1,1000
!print*,'strc',(strc(i1,i2),i2=1,4)
!enddo
!!!!! production of the lone pairs and radical part of the structures starts !!!!!!
if(nlast.ne.0.or.nnlp.ne.0)then
i=0
j5=0
ll=0
do ii=1,alstr
k5=0
i2=0
if(bonds.ne.0)then
do i1=1,nao
do k1=1,bonds*2
if(i1.eq.strct(ii,k1)) goto 530
enddo
i2=i2+1
num(ii,i2)=i1

530 enddo
else
do i1=1,nao
i2=i2+1
num(ii,i2)=i1
enddo
endif
allp=i2

!print*,'num*****',nnlp,(num(ii,i2),i2=1,allp)

!!!!! production of the lone pair part of the structures starts !!!!!!

if(nnlp.ne.0)then
!print*,'cov:lp',nnlp,allp
j=0
m=0
do i1=1,allp
!if(nnlp.eq.3)then
!print*,'allp',i1
!endif
j=1
n(1)=i1
!if(nnlp.eq.3)print*,'sourav1'
if (j.eq.nnlp) then
i=i+1
j4=0
do j1=nnlp,1,-1
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
!print*,(strc(i,l),l=1,nao)
goto 399
endif

do i2=i1+1,allp
!if(nnlp.eq.3)print*,'sourav2',i2
j=2
if(i2.eq.n(1))goto 401
n(2)=i2
if (j.eq.nnlp) then
i=i+1
j4=0
do j1=1,nnlp
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
!print*,(strc(i,l),l=1,nao)
goto 401
endif

do i3=i2+1,allp
!if(nnlp.eq.3)print*,'sourav2',i3
j=3
do k=1,2
if(i3.eq.n(k))goto 402
enddo
n(3)=i3
!print*,'n3',n(3)
if (j.eq.nnlp) then
i=i+1
j4=0
!if(nnlp.eq.3)print*,'sourav3',num(1,n(3)),i,nnlp,ii
!print*,(strc(i,j1),j1=1,2)
do j1=1,nnlp
!print*,'jjj',j
!print*,num(ii,n(j1))
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
!print*,'strc',strc(i,j4)
enddo
enddo

!if(nnlp.eq.3)print*,'sourav4'
if(bonds.ne.0)then
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
endif
!print*,'sourav_strc',(strc(i,l),l=1,j4)
goto 402
endif

do i4=i3+1,allp
j=4
do k=1,3
if(i4.eq.n(k))goto 403
enddo
n(4)=i4
!print*,(n(l),l=1,nao)
if (j.eq.nnlp) then
i=i+1
j4=0
do j1=1,nnlp
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
!print*,'num',i,num(ii,n(j1)),strc(i,j4)
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
!print*,'sourav gadha'
!print*,'strc',(strc(i,l),l=1,nae)
goto 403
endif

!print*,'sourav4'
do i5=i4+1,allp
j=5
do k=1,4
if(i5.eq.n(k))goto 404
enddo
n(5)=i5
if (j.eq.nnlp) then
i=i+1
j4=0
do j1=1,nnlp
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
!print*,(strc(i,l),l=1,nao)
goto 404
endif

do i6=i5+1,allp
j=6
do k=1,5
if(i6.eq.n(k))goto 405
enddo
n(6)=i6
if (j.eq.nnlp) then
i=i+1
j4=0
do j1=1,nnlp
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
!print*,(strc(i,l),l=1,nao)
goto 405
endif

do i7=i6+1,allp
j=7
do k=1,6
if(i7.eq.n(k))goto 406
enddo
n(7)=i7
if (j.eq.nnlp) then
i=i+1
j4=0
do j1=1,nnlp
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 406
endif

do i8=i7+1,allp
j=8
do k=1,7
if(i8.eq.n(k))goto 407
enddo
n(8)=i8
if (j.eq.nnlp) then
i=i+1
j4=0
do j1=1,nnlp
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 407
endif

do i9=i8+1,allp
j=9
do k=1,8
if(i9.eq.n(k))goto 408
enddo
n(9)=i9
if (j.eq.nnlp) then
i=i+1
j4=0
do j1=1,nnlp
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 408
endif

do i10=i9+1,allp
j=10
do k=1,9
if(i10.eq.n(k))goto 409
enddo
n(10)=i10
if (j.eq.nnlp) then
i=i+1
j4=0
do j1=1,nnlp
do j2=1,2
j4=j4+1
strc(i,j4)=num(ii,n(j1))
enddo
enddo
j2=0
do j3=1,bonds*2
j2=j3+j4
strc(i,j2)=strct(ii,j3)
enddo
goto 409
endif

409 enddo
408 enddo
407 enddo
406 enddo
405 enddo
404 enddo
403 enddo
402 enddo
401 enddo
399 enddo
endif
if(bonds.eq.0)j2=j4
!do j1=ll+1,i
!print*,'strccc',j1,(strc(j1,i2),i2=1,j2)
!enddo

if(nlast.ne.0)then
do j1=ll+1,i
l=0
do i3=1,allp
do i2=1,nnlp*2,2
!print*,'strcnum',strc(j1,i2),num(ii,i3)
if(strc(j1,i2).eq.num(ii,i3))goto 439
enddo
l=l+1
num1(l)=num(ii,i3)
439 enddo
!print*,'num1',(num1(i1),i1=1,l)

do i2=1,l
j=1
nn(1)=num1(i2)
if(nlast.eq.j)then
j5=j5+1
l4=0
do l2=1,j2
istrc(j5,l2)=strc(j1,l2)
enddo
do l2=j2+1,nae
l4=l4+1
istrc(j5,l2)=nn(l4)
enddo
!print*,'istrc',(istrc(j5,l2),l2=1,nae)
goto 501
endif

do i3=i2+1,l
j=2
nn(2)=num1(i3)
if(nlast.eq.j)then
j5=j5+1
l4=0
do l2=1,j2
istrc(j5,l2)=strc(j1,l2)
enddo
do l2=j2+1,nae
l4=l4+1
istrc(j5,l2)=nn(l4)
enddo
goto 502
endif

do i4=i3+1,l
j=3
nn(3)=num1(i4)
if(nlast.eq.j)then
j5=j5+1
l4=0
do l2=1,j2
istrc(j5,l2)=strc(j1,l2)
enddo
do l2=j2+1,nae
l4=l4+1
istrc(j5,l2)=nn(l4)
enddo
goto 503
endif

do i5=i4+1,l
j=4
nn(4)=num1(i5)
if(nlast.eq.j)then
j5=j5+1
l4=0
do l2=1,j2
istrc(j5,l2)=strc(j1,l2)
enddo
do l2=j2+1,nae
l4=l4+1
istrc(j5,l2)=nn(l4)
enddo
goto 504
endif


do i6=i5+1,l
j=5
nn(5)=num1(i6)
if(nlast.eq.j)then
j5=j5+1
l4=0
do l2=1,j2
istrc(j5,l2)=strc(j1,l2)
enddo
do l2=j2+1,nae
l4=l4+1
istrc(j5,l2)=nn(l4)
enddo
goto 505
endif


do i7=i6+1,l
j=6
nn(6)=num1(i7)
if(nlast.eq.j)then
j5=j5+1
l4=0
do l2=1,j2
istrc(j5,l2)=strc(j1,l2)
enddo
do l2=j2+1,nae
l4=l4+1
istrc(j5,l2)=nn(l4)
enddo
goto 506
endif


506 enddo
505 enddo
504 enddo
503 enddo
502 enddo
501 enddo


enddo
endif

ll=i
enddo
if(nlast.ne.0)then
alstr=j5

strc(x,y)=0


do i1=1,alstr
do k4=1,nae
strc(i1,k4)=istrc(i1,k4)
enddo
enddo
endif

!print*,'j5j5',j5
if (nnlp.ne.0) alstr=i
if (nlast.ne.0) alstr=j5
!do i1=1,alstr
!print*,(strc(i1,k4),k4=1,nae)
!enddo



endif



990 format(30I3)

!!!!!!! new structures generation part ends !!!!!!!!!!!!!!!!!!!!!!!


write(7,902)'number of lone pairs=',nnlp
write(7,*)'                              '
do i=1,alstr
do i1=1,nae
!print*,'sourav1',i,i1,alstr,nae,nlp
fullionstr(i,i1)=strc(i,i1)+niao
enddo
!print*,'fullionstr',i,(fullionstr(i,i1),i1=1,nae)
enddo
deallocate(strc)

!do j1=1,rns
!print*,j1,')',(allionstr(j1,j2),j2=1,nae)
!enddo
!print*,'***********************************'
d=nao
e=nao-nnlp
f=d-e
elporb=nae-nnlp*2
c=factorial(d)/(factorial(e)*factorial(f))
call wigner(elporb,wig2)

d=0
e=0
f=0
d=nao-nnlp
e=nae-2*nnlp
f=(nao-nnlp)-(nae-2*nnlp)
c1=factorial(d)/(factorial(e)*factorial(f))
perm_nstr=wig2*c1*c
write(7,307)' You have ',alstr,' ionic structures of the &
set of ',nnlp,' lone pair'
write(7,307)'among which',wig2*c1*c,'  structures are permissible mathematically'
write(7,*)'*******************************************************************&
*************************'
307 format(2x,a,I7,a,I3,a)
if(symm.eq.1)then

if(sig_sym_flg.eq.1)call symmetry_cal_sig(nlp,fullionstr,alstr,symsc,symq,nssym)
if(sig_sym_flg.ne.1)call symmetry_cal_pi(nlp,fullionstr,alstr,symsc,symq,nssym)

nnn=0
do ijk=nssym,1,-1
do j2=1,alstr
if(symq(j2).eq.ijk)then
nnn=nnn+1
if(niao.eq.0)then
write(7,900),'ion structure',nnn,')',qulsym(j2),(fullionstr(j2,j1),j1=1,nao-nlp+nlp*2)
endif
if(niao.ne.0)then
write(7,901),'ion structure',nnn,')',qulsym(j2),1,':',niao,(fullionstr(j2,j1),j1=1,nao-nlp+nlp*2)
endif
endif
210 enddo
write(7,*)'******************************************************'
enddo
else
do j2=1,alstr
if(niao.eq.0)then
write(7,903),'ion structure',j2,')',(fullionstr(j2,j1),j1=1,nao-nlp+nlp*2)
endif
if(niao.ne.0)then
write(7,904),'ion structure',j2,')',1,':',niao,(fullionstr(j2,j1),j1=1,nao-nlp+nlp*2)
endif
enddo
write(7,*)'******************************************************'
endif
900 format(a,I5,a,x,I3,3x,30(I5))
901 format(a,I5,a,x,I3,3x,I3,x,a,I3,30(I5))
903 format(a,I5,a,x,30(I5))
904 format(a,I5,a,x,I3,x,a,I3,30(I5))

902 format(30x,a,x,I5)


write(7,*)''
write(7,*)'--------------------------------------------------------------------------'
write(7,*)''

jj=alstr
nos=i
nofs=jjj
is=jj
write(10,*)
write(23,*)
write(10,*)'----------------  Ionic Structures  -------------'
write(23,*)'----------------  Ionic Structures  -------------'
write(10,312),perm_nstr,' ionic structures of ',nnlp,' lone pair'
write(23,312),perm_nstr,' ionic structures of ',nnlp,' lone pair'
write(10,913),'IAB','NNB','SBB','fqual'
write(23,913),'IAB','NNB','SBB','fqual'
913 format(x,a,x,a,x,a,x,a)

if(perm_nstr.eq.is.or.nfset.ne.0) call str_selection(fullionstr,nnlp,alstr,perm_nstr)
if(perm_nstr.ne.is.and.nfset.eq.0)then

230 format(a,20I5)
!do i=1,nstr7
!print*,(str5(i,j),j=1,nae)
!enddo
do i=1,nstr7
ll=0
do j=nlp*2+1,nae-nlast,2
ll=ll+1
str_p(i,ll)=prime_num(str5(i,j))*prime_num(str5(i,j+1))
enddo
if(nlast.ne.0)then
do j=nae-nlast+1,nae
ll=ll+1
str_p(i,ll)=prime_num(str5(i,j))
enddo
endif
enddo
strbn=ll

!do i=1,nstr7
!print*,(str_p(i,j),j=1,strbn)
!enddo

do i=1,jj
ll=0
do j=nnlp*2+1,nae-nlast,2
ll=ll+1
ion_str_p(i,ll)=prime_num(fullionstr(i,j))*prime_num(fullionstr(i,j+1))
space(i)=ion_str_p(i,ll)
enddo
if(nlast.ne.0)then
do j=nae-nlast+1,nae
ll=ll+1
ion_str_p(i,ll)=prime_num(fullionstr(i,j))
space(i)=space(i)*prime_num(fullionstr(i,j))
enddo
endif
istrbn=ll
enddo

!do i=1,jj
!print*,i,(ion_str_p(i,j),j=1,istrbn),space(i)
!enddo

lll=0
llll=0
do i=1,jj
do i1=1,nstr7
ll=0
do i3=1,istrbn
do i2=1,strbn
!print*,ion_str_p(i,i3),str_p(i1,i2)
if(ion_str_p(i,i3).eq.str_p(i1,i2))then
ll=ll+1
goto 721
endif
enddo
721 enddo
if(ll.eq.istrbn)then
lll=lll+1
covris(lll)=i
goto 717
endif
enddo 
717 enddo 

llll=0
do i=1,jj
do i1=1,lll
if(i.eq.covris(i1))goto 718
enddo
llll=llll+1
covnris(llll)=i
!print*,'covnris(llll)',llll,covnris(llll)
718 enddo

ll=0
do i=lll+1,lll+llll
ll=ll+1
covris(i)=covnris(ll)
enddo

do i=1,lll+llll
!print*,'1111111111',i
ll=0
if(space(covris(i)).ne.0)then
!print*,'iiiiisourav',i
do i1=i,lll
l2=0
do i2=1,nnlp*2,2
do i3=1,nnlp*2,2
if(fullionstr(covris(i),i2).eq.fullionstr(covris(i1),i3))l2=l2+1
enddo
enddo

if(space(covris(i)).eq.space(covris(i1)).and.l2.eq.nnlp)then
ll=ll+1
strn(ll)=covris(i1)
do j=1,nae
str3(ll,j)=fullionstr(covris(i1),j)
enddo
!print*,i1,ll,space(covris(i1)),(str3(ll,j),j=1,nae)
if(i.ne.i1)space(covris(i1))=0
endif

enddo
!print*,'wig2',wig2,ll,llll
if(ll.eq.wig2)then
call quality_factor(nnlp,str3,ll,quality_fac,str_quality_1,str_quality_2,bondq)
do m119=1,ll
if(niao.eq.0)then
write(10,800),str_quality_1(m119),bondq(m119),str_quality_2(m119)&
,quality_fac(m119),'|',(str3(m119,m20),m20=1,nae)
endif
if(niao.gt.1)then
write(10,801),str_quality_1(m119),bondq(m119),str_quality_2(m119)&
,quality_fac(m119),'|',1,':',niao,(str3(m119,m20),m20=1,nae)
endif
if(niao.eq.1)then
write(10,809),str_quality_1(m119),bondq(m119),str_quality_2(m119)&
,quality_fac(m119),'|',1,1,(str3(m119,m20),m20=1,nae)
endif
enddo
else
l=0
do i1=1,llll
!print*,'space(covris(i)),covnris(i1)',i1,space(covris(i)),space(covnris(i1))
l2=0
do i2=1,nnlp*2,2
do i3=1,nnlp*2,2
if(fullionstr(covris(i),i2).eq.fullionstr(covnris(i1),i3))l2=l2+1
enddo
enddo
if(space(covris(i)).eq.space(covnris(i1)).and.l2.eq.nnlp)then
l=l+1
ll=ll+1
strn(ll)=covnris(i1)
do j=1,nae
str4(l,j)=fullionstr(covnris(i1),j)
enddo
do j=1,nae
str3(ll,j)=fullionstr(covnris(i1),j)
enddo
!print*,'str3ll',ll,covnris(i1),(str3(ll,j),j=1,nae)
if(i.ne.i1+lll)space(covris(lll+i1))=0
endif
enddo

!do i1=1,ll
!print*,'str3',(str3(i1,j),j=1,nae)
!enddo



if(ll.eq.totistr)then
if(l.gt.1) call quality_factor(nnlp,str4,l,quality_fac,str_quality_1,str_quality_2,bondq)
if(l.gt.1) call qult_str_arrange(nnlp,str4,l,quality_fac,str2,q_fac)
l1=0
do i1=l+1,ll
l1=l1+1
do j=1,nae
str3(i1,j)=str2(l1,j)
enddo
enddo
call quality_factor(nnlp,str3,ll,quality_fac,str_quality_1,str_quality_2,bondq)
call vector_rep(nnlp,str3,ll,fvec)
call write_xmi_new_2(nnlp,wig2,str3,ll,q_fac,quality_fac,fvec) 
endif
endif
endif
enddo

endif

312 format(I3,x,a,x,I3,x,a)
 enddo

tnis=nofs

800 format(I3,2x,I3,2x,I3,2x,I3,x,a,x,25I4)
801 format(I3,2x,I3,2x,I3,2x,I3,x,a,x,I2,a,I2,x,25I4)
809 format(I3,2x,I3,2x,I3,2x,I3,x,a,x,I3,I3,x,25I4)

deallocate(strct)
deallocate(istrc)
deallocate(n)
deallocate(nn)
deallocate(num)
deallocate(num1)

print*,'exit ion_struc'

return
end subroutine ion_struc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function factorial(n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,j,n
real*8::factorial

!print*,'enter factorial'
!print*,'n,n
j=1
do i=n,1,-1
j=j*i
enddo
factorial=j
!print*,'fact',j

n=0
!print*,'exit factorial'
return
end function factorial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function dfactorial(n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,j,n
real*8::dfactorial

!print*,'enter dfactorial'
!print*,'n,n
j=1
do i=n,1,-2
j=j*i
enddo
dfactorial=j
!print*,'fact',j

n=0
print*,'exit dfactorial'
return
end function dfactorial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wigner(nnae,wig2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine calculates the number of permisible structures by wigner's
!!! theorem
use commondat
implicit none

integer::wig2,nm,s,i,j,k,l,b,m,nnae
real*8::a,factorial,ss


!print*,'enter wigner2'
j=0
m=mod(nnae,2)
if(m.eq.1)then
ss=0.5
a=2.0*ss+1.0
do i=1,nnae,2
j=j+1
enddo
k=j+1
l=j-1
wig2=a*factorial(nnae)/(factorial(k)*factorial(l))
endif

if(m.eq.0)then
if(mult.eq.1)s=0
if(mult.eq.3)s=1
b=2*s+1
k=(nnae/2)+1+s
l=(nnae/2)-s

wig2=b*factorial(nnae)/(factorial(k)*factorial(l))
endif


return
end subroutine wigner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine close_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

close(7)

end subroutine close_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_rumer_xmi(nl,str,nstr,rumer,rumer_rad,q_fac)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig

integer::i,i7,nstr,m19,m20,nl,rumer(15000),rumer_rad(15000),q_fac(15000),str_quality_1(15000),bondq(15000),allrum,&
str(15000,20),str_quality_2(15000),tqlty,bqlty,sqlty,tnqs,nssym,qulsym(15000),symq(15000),col(1000),sigsym(15000),tnqs_sig&
,rumstr(500,20)
real*8::ovlp
Double Precision::D(1000)

print*,'enter write_rumer_xmi',nstr
allrum=1
tqlty=0
bqlty=0
sqlty=0
i7=0
do m19=1,nstr
if(rumer(m19)*rumer_rad(m19).eq.1)then
i7=i7+1
col(i7)=m19
do m20=1,nae
rumstr(i7,m20)=str(m19,m20)
enddo
if(nfset.eq.3.or.nfset.eq.5)goto 200
if(niao.eq.0)then
write(10,914),str_quality_1(m19),bondq(m19),str_quality_2(m19),'|',(str(m19,m20),m20=1,nae)
endif
if(niao.gt.1)then
write(10,915),str_quality_1(m19),bondq(m19),str_quality_2(m19),'|',1,':',niao,(str(m19,m20),m20=1,nae)
endif
if(niao.eq.1)then
write(10,916),str_quality_1(m19),bondq(m19),str_quality_2(m19),'|',1,1,(str(m19,m20),m20=1,nae)
endif
tqlty=tqlty+str_quality_1(m19)
bqlty=bqlty+bondq(m19)
sqlty=sqlty+str_quality_2(m19)
endif
200 enddo
open(unit=31,file='Rumer_Sets.dat',status='unknown')
if(allrum.eq.1)call All_Rumer_set(rumstr,i7,nl)
if(nfset.eq.5)stop
if(nfset.eq.3)goto 201
if(ovopt.eq.1)then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(10,912)'Overlap of this set of the structures =',1.0-ovlp
endif
write(10,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty

!write(10,910)'intra bond quality','sym break quality','nnbond quality'
!write(10,911),tqlty,sqlty,bqlty


914 format(I3,x,I3,x,I3,x,a,x,25I4)
915 format(I3,x,I3,x,I3,x,a,x,I1,a,I1,x,25I4)
916 format(I3,x,I3,x,I3,x,a,x,I3,I3,x,25I4)
910 format(a,a,I3,x,a,I3,x,a,I3)
912 format(a,3x,F10.3)
!910 format(a,x,a,x,a,x,a)
!911 format(15x,I3,7x,I3,7x,I3)
print*,'exit write_rumer_xmi'
    
!201 close(21)
201 return
end subroutine write_rumer_xmi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine All_Rumer_set(rumstr,setno,nlonep)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! N! number of Runmer set generates here 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/orb1/orbs1,rstr,nlpr
integer::j,i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,n1,k,nlpr
integer::setno,nlonep,rstr(500,20),permutation(100),orbs1(20),orbs2(20),rumstr(500,20)

nlpr=nlonep
print*,'enter All_Rumer_set',setno
do i=1,setno
do i1=1,nae
rstr(i,i1)=rumstr(i,i1)
enddo
enddo
i=0
do j=nlonep*2+1,nae
i=i+1
orbs2(i)=rumstr(1,j)
enddo
n1=i

!print*,'orns_rum',(orbs2(i),i=1,n1)

i=100
do j=1,n1
if(orbs2(j).lt.i)i=orbs2(j)
enddo

orbs1(1)=i
k=1
do j=1,20
i=i+1
!print*,i
do i1=1,n1
if(i.eq.orbs2(i1))then
k=k+1
orbs1(k)=i
print*,i
endif
enddo
enddo
!print*,(orbs1(i),i=1,n1)

j=0
do i1=1,n1
permutation(1)=orbs1(i1)

do i2=1,n1
if(permutation(1).eq.orbs1(i2))goto 201
permutation(2)=orbs1(i2)
if(n1.eq.2)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 201
endif

do i3=1,n1
do i=1,2
if(permutation(i).eq.orbs1(i3))goto 202
enddo
permutation(3)=orbs1(i3)
if(n1.eq.3)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 202
endif

do i4=1,n1
do i=1,3
if(permutation(i).eq.orbs1(i4))goto 203
enddo
permutation(4)=orbs1(i4)
if(n1.eq.4)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 203
endif

do i5=1,n1
do i=1,4
if(permutation(i).eq.orbs1(i5))goto 204
enddo
permutation(5)=orbs1(i5)
if(n1.eq.5)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 204
endif

do i6=1,n1
do i=1,5
if(permutation(i).eq.orbs1(i6))goto 205
enddo
permutation(6)=orbs1(i6)
if(n1.eq.6)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 205
endif
do i7=1,n1
do i=1,6
if(permutation(i).eq.orbs1(i7))goto 206
enddo
permutation(7)=orbs1(i7)
if(n1.eq.7)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 206
endif

do i8=1,n1
do i=1,7
if(permutation(i).eq.orbs1(i8))goto 207
enddo
permutation(8)=orbs1(i8)
if(n1.eq.8)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 207
endif

do i9=1,n1
do i=1,8
if(permutation(i).eq.orbs1(i9))goto 208
enddo
permutation(9)=orbs1(i9)
if(n1.eq.9)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 208
endif

do i10=1,n1
do i=1,9
if(permutation(i).eq.orbs1(i10))goto 209
enddo
permutation(10)=orbs1(i10)
if(n1.eq.10)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 209
endif

do i11=1,n1
do i=1,10
if(permutation(i).eq.orbs1(i11))goto 210
enddo
permutation(11)=orbs1(i11)
if(n1.eq.11)then
j=j+1
call rumer(permutation,n1,j,setno,nlpr)
goto 210
endif
210 enddo
209 enddo
208 enddo
207 enddo
206 enddo
205 enddo
204 enddo
203 enddo
202 enddo
201 enddo
200 enddo
print*,'exit All_Rumer_set'

return
end subroutine All_Rumer_set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rumer(permutation,n,j,setno,nl)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! All the Rumer sets are written in the file "Rumer_Sets.dat" varified
!! with the subroutine "Rumer_set_id" and then written in the file
!! "Rumer_Sets_all.dat" with full format of the output.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
common/orb1/orbs1,rstr,nlonep

integer::permutation(100),n,i,i1,i2,j,orbs1(20),rstr(500,20),setno,m20,m19,rstr1(2000,20),nlonep
Double Precision::bond_rum(10000,100),bond_rum1(10000,100),k1,k2
integer::l,ll,jj,nl,rum_count(10000),set_num(100),Rid

!totrum=j
Rumwrite=0
do i=1,setno
do i1=1,nlonep*2
rstr1(i,i1)=rstr(i,i1)
enddo
enddo

print*,'enter rumer'
do i=1,setno
do i1=1,nae
do i2=1,n
if(rstr(i,i1).eq.orbs1(i2))then
!write(21,*),i2,rstr(i,i1),orbs1(i2),permutation(i2)
rstr1(i,i1)=permutation(i2)
!write(21,*),i1,rstr(i,i1)
endif
enddo
enddo
enddo


if(j.eq.1)then
jj=0
goto 103
endif
call Rumer_set_id(rstr1,setno,nlonep,Rid,set_num) 
if(Rid.eq.1)goto 100



103 jj=jj+1

write(31,*)'set number',jj,(permutation(i),i=1,n),setno
write(23,*)'permutation',jj,'>',(permutation(i),i=1,n)
write(23,*)'*****************************************************'
write(23,*)
do m19=1,setno
write(31,914),(rstr1(m19,m20),m20=1,nae)
if(nfset.eq.5)write(10,915),1,':',niao,(rstr1(m19,m20),m20=1,nae)
!if(niao.eq.0)then
!write(23,914),(rstr1(m19,m20),m20=1,nae)
!endif
if(niao.gt.1)then
!write(23,915),1,':',niao,(rstr1(m19,m20),m20=1,nae)
!write(10,915),q_fac(m19),1,':',niao,(str(m19,m20),m20=1,nae)
endif
!if(niao.eq.1)then
!write(23,916),1,1,(rstr1(m19,m20),m20=1,nae)
!endif
enddo
!write(23,*)
if(nfset.eq.5)write(10,*)'set number',jj

900 format(a,I3,10I5)

914 format(x,25I4)
915 format(x,I1,a,I3,x,25I4)
916 format(x,I3,I3,x,25I4)
print*,'exit rumer'
totrum=jj
print*,'totrum',totrum
!stop
100 return
end subroutine rumer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rumer_structures(nl,str,nstr,rumer,rumer_rad)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::kk,k2,k3,k4,k5,k6,k7,k8,k9,nl,nstr,rum,rum1,rum2,rum3,rum4,rrad(100),rrad1,x,y
integer::str(15000,20),rumer(15000),rumer_rad(15000),rumer1(15000,100),rumer2(15000,100)

print*,'enter rumer_structures'

do k6=1,15000
rumer(k6)=0
rumer_rad(k6)=0
enddo
do k3=1,15000
do k6=1,100
rumer1(k3,k6)=0
enddo
enddo
do k3=1,nstr
do k6=1,nae-nl*2-nlast
do k7=nl*2+1,nae-nlast
do k8=1,k6
if(rumer1(k3,k8).eq.str(k3,k7))goto 389
enddo
if(str(k3,k7).gt.rumer1(k3,k6))then
rumer1(k3,k6)=str(k3,k7)
!print*,'**',rumer1(k3,k6)
endif
389 enddo
enddo
enddo
do k3=1,15000
do k6=1,100
rumer2(k3,k6)=0
enddo
enddo

231 format(20I3)
do k3=1,nstr
k7=0
do k6=nae-nl*2-nlast,1,-1
k7=k7+1
rumer2(k3,k7)=rumer1(k3,k6)

!print*,rumer2(k3,k7)
enddo
enddo
do k3=1,nstr
print*,'**',(rumer2(k3,k6),k6=1,nae-nl*2-nlast)
enddo
do k6=1,15000
rumer(k6)=0
enddo
!print*,'k7',k7
do k9=1,nstr
!write(*,231),(str(k9,k6),k6=1,nae)
rum=0
do k6=nl*2+1,nae-nlast,2
do k2=1,k7
k8=k6
if(str(k9,k8).eq.rumer2(k9,k2))rum1=k2
enddo
do k2=1,k7
k8=k6+1
if(str(k9,k8).eq.rumer2(k9,k2))rum2=k2
enddo

do k5=k6+2,nae-nlast,2
do k2=1,k7
k8=k5
if(str(k9,k8).eq.rumer2(k9,k2))rum3=k2
enddo
do k2=1,k7
k8=k5+1
if(str(k9,k8).eq.rumer2(k9,k2))rum4=k2
enddo

if((rum1-rum3)*(rum1-rum4)*(rum2-rum3)*(rum2-rum4).gt.0)then

rum=rum+1
endif
enddo
enddo
kk=0
do k6=1,((nae-nl*2-nlast)/2)-1
kk=kk+k6
enddo

!print*,'rum1,rum2,rum',rum1,rum2,rum
if(rum.eq.kk)then
rumer(k9)=1
!print*,'rumer(k9)',rumer(k9)
endif
enddo

do k6=1,15000
rumer_rad(k6)=1
enddo
if(nlast.ne.0)then
do k3=1,15000
do k6=1,100
rumer2(k3,k6)=0
rumer1(k3,k6)=0
enddo
enddo
do k3=1,nstr
do k6=1,nae-nl*2
do k7=nl*2+1,nae
do k8=1,k6
if(rumer1(k3,k8).eq.str(k3,k7))goto 390
enddo
if(str(k3,k7).gt.rumer1(k3,k6))then
rumer1(k3,k6)=str(k3,k7)
!print*,rumer1(k6)
endif
390 enddo
enddo
enddo
do k3=1,15000
do k6=1,100
rumer2(k3,k6)=0
enddo
enddo
do k3=1,15000
k7=0
do k6=nae-nl*2,1,-1
k7=k7+1
rumer2(k3,k7)=rumer1(k3,k6)

enddo
enddo

do k9=1,nstr
rum=0
do k6=nae-nlast+1,nae
do k2=1,k7
!print*,str(k9,k6),rumer2(k9,k2)
if(str(k9,k6).eq.rumer2(k9,k2))rrad(k6)=k2
enddo
enddo
do k6=nl*2+1,nae-nlast,2
do k2=1,k7
k8=k6
if(str(k9,k8).eq.rumer2(k9,k2))rum1=k2
enddo
do k2=1,k7
k8=k6+1
if(str(k9,k8).eq.rumer2(k9,k2))rum2=k2
enddo
do k4=nae-nlast+1,nae
!print*,'rrad,rum1,rum2',rrad(k4),rum1,rum2
if((rrad(k4)-rum1)*(rrad(k4)-rum2).gt.0)then
rum=rum+1
endif
enddo
enddo
if(rum.lt.((nae-nlast-nl*2)/2)*nlast)then
rumer_rad(k9)=0
!print*,'rumer_rad(k9)',rumer_rad(k9)
endif
enddo


endif

print*,'exit rumer_structures'
return
end subroutine rumer_structures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sym_check(nl,str2,n,sym_str_sl,sym_str_num,numsymset)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig

logical :: fileexists
character(len=35)::inputfilename
integer::i,i1,i2,i3,i4,i5,i6,i7,i8,l,l1,m,bnd,n,nl,numsymset,nsymset,nsymstr,str6(100,20),str2(15000,20)
integer::str_quality_1(15000),str_quality_2(15000),bondq(15000),tqlty,bqlty,sqlty,tnqs,nssym,sym_check_flg,&
qulsym(15000),symq(15000),sigsym(15000),tnqs_sig,sym_str_sl1(1000,100),sym_str_sl(1000,100),sym_str_num(1000),&
sym_str_num1(1000),val,stc

print*, 'enter sym_check'
inputfilename='sym_str_set.dat'
INQUIRE(FILE=TRIM(inputfilename),EXIST=fileexists)
IF (fileexists) THEN
open(unit=43,file='sym_str_set.dat',status='unknown')
ELSE
PRINT*,'SORRY sym_str_set.dat file does not exist'
print*,'please write the file: at the top need to mention the number of sets and then number of structures in each sets'
print*,'please leave a empty line after each segment'
stop
ENDIF



read(43,*)nsymset,nsymstr
stc=0
bnd=(nae-nl*2-nlast)/2

101 do i=1,nsymset
read(43,*)
do i2=1,numsymset
do i4=1,sym_str_num(i2)
sym_str_sl1(i2,i4)=sym_str_sl(i2,i4)
enddo
enddo
do i1=1,nsymstr
read(43,*)(str6(i1,i2),i2=1,nae)
do i3=1,n
m=0
do i4=nl*2+1,nl*2+bnd*2,2
do i7=nl*2+1,nl*2+bnd*2,2
l=0
do i5=i4,i4+1
do i8=i7,i7+1
if(str6(i1,i5).eq.str2(i3,i8))then
l=l+1
goto 109
endif
enddo
109 enddo
if(l.eq.2)then
m=m+1
goto 111
endif
enddo
goto 112
111 enddo
if(m.eq.bnd)then
do i2=1,numsymset
do i4=1,sym_str_num(i2)
if(sym_str_sl1(i2,i4).eq.i3)then
sym_str_sl1(i2,i4)=0
endif
enddo
enddo
endif
112 enddo
enddo

do i2=1,1000
sym_str_num1(i2)=0
enddo

do i2=1,numsymset
do i4=1,sym_str_num(i2)
if(sym_str_sl1(i2,i4).ne.0)val=1
if(sym_str_sl1(i2,i4).eq.0)val=0
sym_str_num1(i2)=sym_str_num1(i2)+val
enddo
enddo

!do i2=1,numsymset
!print*,sym_str_num1(i2),sym_str_num(i2)
!enddo

!do i3=1,numsymset
!do i2=1,sym_str_num(i)
!write(*,200),i,')',(str2(sym_str_sl(i3,i2),i4),i4=1,nae)
!enddo
!enddo

sym_check_flg=1
do i2=1,numsymset
if(sym_str_num(i2)-sym_str_num1(i2).ne.0.and.sym_str_num(i2)-sym_str_num1(i2).ne.sym_str_num(i2))then
sym_check_flg=0
goto 114
endif
enddo
114 do i3=1,nsymstr
write(10,202)(str6(i3,i2),i2=1,nae)
enddo
if(sym_check_flg.eq.1)then
stc=stc+1
write(10,201)'set',i,'is symmetric'
endif
if(sym_check_flg.eq.0)write(10,201)'set',i,'is not symmetric'
write(10,*)

enddo
write(10,*),'Total number of symmetric sets =',stc

200 format (I2,a,2x,20I3)
201 format (a,I5,2x,a)
202 format (20I3)

!do i2=1,numsymset
!print*,(sym_str_sl(i2,i3),i3=1,sym_str_num(i2))
!enddo
!print*,'******************'
!do i2=1,numsymset
!do i3=1,sym_str_num(i2)
!print*,i2,sym_str_sl(i2,i3),(str2((sym_str_sl(i2,i3)),i4),i4=1,nae)
!enddo
!enddo

print*, 'exit sym_check'
stop
end subroutine sym_check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine qult_str_arrange(nl,str2,n,q_fac,str1,q_fac1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig

integer::i1,i5,i6,i7,l,l1,k,m,bnd,nsymstr,ll,nstsymset,sym_str_num_final(1000),sym_str_qual(1000)
integer::sym_str_sl(1000,100),sym_str_sl1(1000),sym_str_sl_final(1000,100),p,nstrt1,&
sym_str_num(1000),sym_str_num1(1000),val,mq_fac,msymq,mloopsymsc,ind,sl,sll
integer::i,i2,i3,i4,i9,m16,m18,m19,q_fac(15000),q_fac1(15000),k6,k7,k8,k9,st_ct(15000),st_ct1(15000),ii,kk
integer::str1(15000,20),str2(15000,20),n,jj,j,j1,nqul1,nl,iii,x,y,z,iij,i8,q_cnt(1000)&
,str_quality_1(15000),str_quality_2(15000),bondq(15000),sigsym(15000),tnqs_sig&
,tqlty,bqlty,sqlty,tnqs,nssym,qulsym(15000),symq(15000),symq_final(15000),ij,ijk,qt,qtg,q&
,mat(100,100),cnt(15000),qt_max,q_max,lll,sym_str_qual1(1000,100),str_cnt1(1000000),nssym1,mm
integer,dimension(:,:),allocatable::str3
integer,dimension(:),allocatable::qual4,qual3,str_cnt,q_fac3,q_fac4,q_fac5,q_fac6,q_fac7,q_fac2,bdq,bdq1,qul1
integer,dimension(:),allocatable::stq1,stq2,stq3,stq4,sl_num,sl_qfac,sl_qfac1

!set_order desides if the symmetric sets arranged according to quality or size
! set_order=0 means quality
! set_order=1 means  smaller to larger sets
! set_order=2 means larger to smaller sets
print*,'set_order',set_order
!stop
print*,'enter qult_str_arrange'
x=50000
y=20
z=1000

allocate(q_fac2(x))
allocate(q_fac3(x))
allocate(q_fac4(x))
allocate(q_fac5(x))
allocate(q_fac6(x))
allocate(q_fac7(x))
allocate(qul1(x))
allocate(qual3(x))
allocate(qual4(x))
allocate(bdq(x))
allocate(bdq1(x))
allocate(str_cnt(z))
allocate(stq1(x))
allocate(stq2(x))
allocate(stq3(x))
allocate(stq4(x))
allocate(sl_num(x))
allocate(sl_qfac(x))
allocate(sl_qfac1(x))

q_fac3(x)=0
q_fac4(x)=0
q_fac5(x)=0
q_fac6(x)=0
q_fac7(x)=0
q_fac2(x)=0
qul1(x)=0
qual3(x)=0
qual4(x)=0
bdq(x)=0
bdq1(x)=0
str_cnt(z)=0
stq1(x)=0
stq2(x)=0
stq3(x)=0
stq4(x)=0
sl_num(x)=0
sl_qfac(x)=0
sl_qfac1(x)=0

do i=1,15000
st_ct1(i)=0
enddo

print*,'nstrtIIIIIIIIIIIIIIIII',nstrt

jj=1
do m19=1,n
write(*,231),'str2:qarng',m19,(str2(m19,m18),m18=1,nae),q_fac(m19),symq(m19)
!print*,'str2',m19,q_fac(m19)
if(m19.eq.1)qul1(1)=q_fac(1)
j=jj
do i=1,j
if(qul1(i).eq.q_fac(m19))goto 373
enddo
jj=jj+1
qul1(jj)=q_fac(m19)
!print*,'qul1',qul1(jj)
373 enddo
nqul1=jj

print*,'nqul1',nqul1


do k6=1,15000
qual3(k6)=0
enddo

print*,'check1'
do k6=1,nqul1
do k7=1,nqul1
do k8=1,k6
if(qual3(k8).eq.qul1(k7))goto 389
enddo
if(qul1(k7).gt.qual3(k6))then
qual3(k6)=qul1(k7)
print*,'**',qual3(k6)
endif
389 enddo
enddo
print*,'check2'

k7=0
do k6=jj,1,-1
k7=k7+1
qual4(k7)=qual3(k6)
enddo

print*,'check3'

do m18=1,15000
do m19=1,nae
str1(m18,m19)=0
enddo
enddo

print*,n,noqult
!!!!! for the non-symmetric system .. if user dont want to take the symmetry in the calculations !!!!!!!!! 
if(symm.eq.0)then
i4=0
do m18=1,nqul1
do m19=1,n
!print*,q_fac(m19),qult(m18)
if(q_fac(m19).eq.qual4(m18))then 
i4=i4+1
do i3=1,nae
str1(i4,i3)=str2(m19,i3)
enddo
q_fac1(i4)=q_fac(m19)
stq1(i4)=str_quality_1(m19)
stq2(i4)=str_quality_2(m19)
bdq(i4)=bondq(m19)
!lfst1(i4)=lfst(m19)
endif
enddo
enddo


endif

!do m18=1,n
!print*,'sourav1',m18,(str2(m18,m19),m19=1,nae),qulsym(m18),symq(m18),q_fac(m18),tnqs,nssym,nqul1
!enddo

!!!!! for the symmetric system .. if user opted to take the symmetry calculations !!!!!!!!! 
if(symm.ne.0)then
do i=1,n
print*,'i,qfac',i,q_fac(i),symq(i)
enddo

print*,'check3'
ind=0
do i=1,n
if(ind.lt.q_fac(i))ind=q_fac(i)
enddo
mq_fac=ind

ind=0
do i=1,n
if(ind.lt.symq(i))ind=symq(i)
enddo
msymq=ind

ind=0
do i=1,n
if(ind.lt.loopsymsc(i))ind=loopsymsc(i)
enddo
mloopsymsc=ind
print*,'mq_fac,msymq,mloopsymsc',mq_fac,msymq,mloopsymsc

bnd=(nae-nl*2-nlast)/2
m=0
ll=0
mm=0
301 do i=1,mq_fac
if(mloopsymsc.ne.0)then
do k=1,mloopsymsc
do j=1,msymq
lll=1
l=0
ll=ll+1
do i1=1,n
do i2=1,mm
if(str_cnt1(i2).eq.i1)goto 307
enddo
if(q_fac(i1).eq.i)then
if(loopsymsc(i1).eq.k)then
if(symq(i1).eq.j)then
lll=0
m=m+1
mm=mm+1
print*,'mmmmm',m,i1,i,j,ll,loopsymsc(i1)
l=l+1

str_cnt1(mm)=i1
sym_str_sl(ll,l)=i1
sym_str_qual(ll)=q_fac(i1)
sym_str_qual1(ll,l)=q_fac(i1)
endif
endif
endif
307 enddo
if(lll.eq.0)sym_str_num(ll)=l

if(lll.eq.1)ll=ll-1
if(m.eq.n)then

nssym1=ll
goto 101
endif
enddo
enddo
endif
if(mloopsymsc.eq.0)then
do j=1,msymq
lll=1
l=0
ll=ll+1
do i1=1,n
do i2=1,mm
if(str_cnt1(i2).eq.i1)goto 407
enddo
if(q_fac(i1).eq.i)then
if(symq(i1).eq.j)then
lll=0
m=m+1
mm=mm+1
print*,'mmmmm',m,i1,i,j,ll,loopsymsc(i1)
l=l+1

str_cnt1(mm)=i1
sym_str_sl(ll,l)=i1
sym_str_qual(ll)=q_fac(i1)
sym_str_qual1(ll,l)=q_fac(i1)
endif
endif
407 enddo
if(lll.eq.0)sym_str_num(ll)=l

if(lll.eq.1)ll=ll-1
if(m.eq.n)then

nssym1=ll
goto 101
endif
enddo
endif
enddo
print*,'check4',m,n
goto 301

101 do i2=1,nssym1
print*,(sym_str_sl(i2,i3),i3=1,sym_str_num(i2)),'|',(sym_str_qual1(i2,i3),i3=1,sym_str_num(i2)),'|',sym_str_qual(i2)
enddo

!print*,'sourav_is_here'
!do i2=1,nssym1
!print*,'sym_str_num',sym_str_num(i2)
!enddo


if(symtype.eq.'check')then
call sym_check(nl,str2,n,sym_str_sl,sym_str_num,nssym1)
endif

!nstsymset=0
!ll=0
!if(nstrt.ne.0)then
!do i1=1,nstrt
!print*,(strt_struc(i1,i5),i5=1,nae)
!do i3=1,n
!!print*,'i3333',i3
!m=0
!do i4=nl*2+1,nl*2+bnd*2,2
!do i7=nl*2+1,nl*2+bnd*2,2
!l=0
!do i5=i4,i4+1
!do i8=i7,i7+1
!!print*,'str6,str2',str6(i1,i5),str2(i3,i8)
!if(strt_struc(i1,i5).eq.str2(i3,i8))then
!l=l+1
!goto 109
!endif
!enddo
!109 enddo
!!print*,'llll',l
!if(l.eq.2)then
!m=m+1
!goto 111
!endif
!enddo
!goto 112
!111 enddo
!!print*,'mmmm',m,bnd
!if(m.eq.bnd)then
!!print*,'i3i3i3',i3
!do i2=1,nssym1
!do i4=1,sym_str_num(i2)
!if(sym_str_sl(i2,i4).eq.i3)then
!ll=ll+1
!sym_str_sl1(ll)=i2
!print*,'sym_str_sl1(ll)',sym_str_sl1(ll)
!endif
!enddo
!enddo
!endif
!112 enddo
!enddo
!!!!! nstsymset=number of starting symmetry set
!nstsymset=ll

!print*,'ll',ll,sym_str_sl1(1)
!
!
!do i2=1,nssym1
!write(*,*),'sym_str_sl',(sym_str_sl(i2,i3),i3=1,sym_str_num(i2))
!enddo
!
!do i4=1,nstsymset
!do i3=1,sym_str_num(sym_str_sl1(i4))
!sym_str_sl_final(i4,i3)=sym_str_sl(sym_str_sl1(i4),i3)
!enddo
!sym_str_num_final(i4)=sym_str_num(sym_str_sl1(i4))
!enddo 

!i9=0
!do i4=1,nstsymset
!do i3=1,sym_str_num(sym_str_sl1(i4))
!i9=i9+1
!str_cnt(i9)=sym_str_sl(sym_str_sl1(i4),i3)
!symq_final(i9)=symq(sym_str_sl(sym_str_sl1(i4),i3))
!do j=1,nae
!str1(i9,j)=str2(sym_str_sl(sym_str_sl1(i4),i3),j)
!enddo
!q_fac1(i9)=q_fac(sym_str_sl(sym_str_sl1(i4),i3))
!stq1(i9)=str_quality_1(sym_str_sl(sym_str_sl1(i4),i3))
!stq2(i9)=str_quality_2(sym_str_sl(sym_str_sl1(i4),i3))
!bdq(i9)=bondq(sym_str_sl(sym_str_sl1(i4),i3))
!enddo
!enddo



!endif

l=0
do i=1,1000
do i2=1,nssym1
!if(nstrt.ne.0)then
!do i4=1,nstsymset
!if(i2.eq.sym_str_sl1(i4))goto 115
!enddo
!endif
if(sym_str_qual(i2).eq.i)then
l=l+1
sym_str_num_final(l)=sym_str_num(i2)
!print*,'sym_str_num_final(l)',sym_str_num_final(l)
do i3=1,sym_str_num(i2)
sym_str_sl_final(l,i3)=sym_str_sl(i2,i3)
!print*,'sym_str_sl_final(l,i3)',sym_str_sl_final(l,i3)
enddo
endif
115 enddo
enddo
!endif

!do i2=1,nssym1
!print*,'sym_str_num_final(i2)',sym_str_num_final(i2)
!enddo

if(set_order.eq.1.or.set_order.eq.2)then
do l=1,nssym1
sym_str_num(l)=sym_str_num_final(l)
do i3=1,sym_str_num(l)
sym_str_sl(l,i3)=sym_str_sl_final(l,i3)
enddo
enddo

endif


!!!bigger set is going top
if(set_order.eq.2)then
l=ll
do i=100,1,-1
do i2=1,nssym1
!if(nstrt.ne.0)then
!do i4=1,nstsymset
!if(i2.eq.sym_str_sl1(i4))goto 113
!enddo
!endif
if(sym_str_num(i2).eq.i)then
l=l+1
sym_str_num_final(l)=sym_str_num(i2)
do i3=1,sym_str_num(i2)
sym_str_sl_final(l,i3)=sym_str_sl(i2,i3)
enddo
endif
113 enddo
enddo
endif


!!!smaller sets is going top
if(set_order.eq.1)then

l=ll
do i=1,100
do i2=1,nssym1
!if(nstrt.ne.0)then
!do i4=1,nstsymset
!if(i2.eq.sym_str_sl1(i4))goto 114
!enddo
!endif
if(sym_str_num(i2).eq.i)then
l=l+1
sym_str_num_final(l)=sym_str_num(i2)
do i3=1,sym_str_num(i2)
sym_str_sl_final(l,i3)=sym_str_sl(i2,i3)
enddo
endif
114 enddo
enddo

endif


print*,'*********************************'
print*,'nssym1',nssym1
do i2=1,nssym1
print*,'i2',i2
print*,'sym_str_num_final(i2)',sym_str_num_final(i2)
write(*,906),(sym_str_sl_final(i2,i3),i3=1,sym_str_num_final(i2))
enddo

906 format (50I3)

i9=0
do i1=1,nssym1
do i2=1,sym_str_num_final(i1)
i9=i9+1
str_cnt(i9)=sym_str_sl_final(i1,i2)
do j=1,nae
str1(i9,j)=str2(sym_str_sl_final(i1,i2),j)
enddo
q_fac1(i9)=i1
!q_fac1(i9)=q_fac(sym_str_sl_final(i1,i2))
stq1(i9)=str_quality_1(sym_str_sl_final(i1,i2))
stq2(i9)=str_quality_2(sym_str_sl_final(i1,i2))
bdq(i9)=bondq(sym_str_sl_final(i1,i2))
enddo
enddo

do i9=1,n
write(*,231),'str2str2',i9,(str1(i9,m18),m18=1,nae),q_fac1(i9),bdq(i9)
enddo
231 format(a,30I3)
endif
do i=1,n
str_quality_1(i)=stq1(i)
str_quality_2(i)=stq2(i)
bondq(i)=bdq(i)
enddo

do i=1,n
write(*,231),'symmstr',i,(str1(i,m18),m18=1,nae),q_fac1(i),bdq(i),str_quality_1(i),str_quality_2(i)
enddo

!!! if user want some symmetric set must be availeble stats bellow !!!

allocate(str3(x,y))
str3(x,y)=0
if(nstrt.ne.0.and.symm.eq.1)then
do j=1,n
do k=1,nstrt
if(nl.ne.0)then
sl=0
do i=1,nl*2,2
do i1=1,nl*2,2
if(str1(j,i).eq.strt_struc(k,i1))then
sl=sl+1
endif
enddo
enddo
if (sl.eq.nl)goto 372
goto 471
endif

372 sll=0
do i=nl*2+1,nae-nlast,2
do i1=nl*2+1,nae-nlast,2
sl=0
do i3=i,i+1
do i4=i1,i1+1
if(str1(j,i3).eq.strt_struc(k,i4))then
sl=sl+1
if(sl.eq.2)then
sll=sll+1
goto 473
endif
endif
enddo
enddo
enddo
473 enddo
if(sll.eq.(nae-nlast-nl*2)/2)goto 374
goto 471

374 if(nlast.ne.0)then
sl=0
do i=nae-(nlast-1),nae
do i1=nae-(nlast-1),nae
if(str1(j,i).eq.strt_struc(k,i1))then
sl=sl+1
endif
enddo
enddo
if(sl.eq.nlast)goto 375
goto 471
endif

375 sl_qfac(k)=q_fac1(j)
print*,'sl_qfac',sl_qfac(k)
471 enddo
enddo

l=1
sl_qfac1(l)=sl_qfac(l)
do i=1,nstrt
p=sl_qfac(i)
if(sl_qfac(i-1).eq.p.or.i.eq.1) goto 477
l=l+1
sl_qfac1(l)=p
print*,'sl_qfac(i)',sl_qfac(i)
477 enddo
nstrt1=l

print*,'lllll',l
do i=1,nstrt1
print*,'sl_qfac1(l)',sl_qfac1(i)
enddo

j1=0
l=0
do j=1,n
do i=1,nstrt1
if(sl_qfac1(i).eq.q_fac1(j))then
l=l+1
sl_num(l)=j
print*,'sl_num(l)',sl_num(l)
goto 376
endif
enddo
!j1=j1+1
!print*,'jjjjjjjjj',j,j1
!do i1=1,nae
!str3(j1,i1)=str1(j,i1)
!enddo
!q_fac3(j1)=q_fac1(j)
!stq1(j1)=str_quality_1(j)
!stq2(j1)=str_quality_2(j)
!bdq(j1)=bondq(j)
376 enddo

do i=1,j1
print*,'i,str3(i,i1)',i,(str3(i,i1),i1=1,nae)
enddo
do j=1,n
do i=1,l
if(sl_num(i).ne.j)goto 377
j1=j1+1
do i1=1,nae
str3(j1,i1)=str1(j,i1)
enddo
q_fac3(j1)=q_fac1(j)
stq1(j1)=str_quality_1(j)
stq2(j1)=str_quality_2(j)
bdq(j1)=bondq(j)
377 enddo
enddo

do j=1,n
do i=1,l
if(sl_num(i).eq.j)goto 378
enddo
j1=j1+1
do i1=1,nae
str3(j1,i1)=str1(j,i1)
enddo
q_fac3(j1)=q_fac1(j)
stq1(j1)=str_quality_1(j)
stq2(j1)=str_quality_2(j)
bdq(j1)=bondq(j)
378 enddo


do j=1,n
print*,'j,rearranged',j,(str3(j,i1),i1=1,nae),q_fac3(j)
enddo
do i=1,n
do i2=1,nae
str1(i,i2)=str3(i,i2)
enddo
q_fac1(i)=q_fac3(i)
str_quality_1(i)=stq1(i)
str_quality_2(i)=stq2(i)
bondq(i)=bdq(i)
enddo
endif
deallocate(str3)
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! below part works when user wish to have some structures always in the top of the list !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(str3(x,y))
str3(x,y)=0

if(nstrt.ne.0.and.symm.eq.0)then
jj=0
if(nl.ne.0)then
jj=0
do i3=1,nstrt
ii=0
do i=1,nl*2,2 
do i2=1,nl*2,2
if(strt_struc(i3,i).eq.str1(1,i2))then
ii=ii+1
endif
enddo
enddo
if(ii.eq.nl)then
jj=jj+1
st_ct(jj)=i3
endif
enddo
endif

print*,'***********nl,jj***************',nl,jj
if(jj.ne.0)then
kk=0
do i=1,n
!print*,'str1',(str1(i,m19),m19=1,nae)
do k8=1,jj
iii=0
k9=st_ct(k8)
!print*,'strt_struc',(strt_struc(k9,m19),m19=1,nae)
do i2=nl*2+1,nae-nlast,2
do k6=nl*2+1,nae-nlast,2
ii=0
do i4=i2,i2+1
do i3=k6,k6+1
!print*,'k9,i3,i4,strt_struc(k9,i3),str1(n,i4)',k9,i3,i4,strt_struc(k9,i3),str1(n,i4)
if(strt_struc(k9,i3).eq.str1(i,i4))then
ii=ii+1
!print*,'ii',ii
endif
enddo
enddo
if(ii.eq.2)then
iii=iii+ii
!print*,'ii,iii**',ii,iii,nae-nl*2-nlast
if(iii.eq.nae-nl*2-nlast)then
kk=kk+1
st_ct1(kk)=i
!print*,'st_ct1(kk)',st_ct1(kk)
if(kk.eq.jj)goto 340
goto 360
endif
goto 350
endif
enddo
350 enddo
enddo
360 enddo


340 jj=0
do i=1,kk
jj=jj+1
do m18=1,nae
str3(jj,m18)=str1(st_ct1(i),m18)
enddo
q_fac2(jj)=q_fac1(st_ct1(i))
stq3(jj)=stq1(st_ct1(i))
stq4(jj)=stq2(st_ct1(i))
bdq1(jj)=bdq(st_ct1(i))
enddo

do i=1,n
do i2=1,kk
if(i.eq.st_ct1(i2))goto 370
enddo
jj=jj+1
do m18=1,nae
str3(jj,m18)=str1(i,m18)
enddo
q_fac2(jj)=q_fac1(i)
stq3(jj)=stq1(i)
stq4(jj)=stq2(i)
bdq1(jj)=bdq(i)
370 enddo

print*,'souravooooo'
do i=1,n
q_fac1(i)=0
str_quality_1(i)=0
str_quality_2(i)=0
bondq(i)=0
do i2=1,nae
str1(i,i2)=0
enddo
enddo

do i=1,n
q_fac1(i)=q_fac2(i)
str_quality_1(i)=stq3(i)
str_quality_2(i)=stq4(i)
bondq(i)=bdq1(i)
do i2=1,nae
str1(i,i2)=str3(i,i2)
enddo
enddo

endif


if(nl.eq.0)then
kk=0
do i=1,n
!print*,'str1',(str1(i,m19),m19=1,nae)
do k9=1,nstrt
iii=0
print*,'strt_struc',(strt_struc(k9,m19),m19=1,nae)
do i2=1,nae-nlast,2
do k6=1,nae-nlast,2
ii=0
do i4=i2,i2+1
do i3=k6,k6+1
!print*,'k9,i3,i4,strt_struc(k9,i3),str1(n,i4)',k9,i3,i4,strt_struc(k9,i3),str1(n,i4)
if(strt_struc(k9,i3).eq.str1(i,i4))then
ii=ii+1
!print*,'ii',ii
endif
enddo
enddo
if(ii.eq.2)then
iii=iii+ii
!print*,'ii,iii**',ii,iii,nae-nl*2-nlast
if(iii.eq.nae-nlast)then
kk=kk+1
st_ct1(kk)=i
!print*,'st_ct1(kk)',st_ct1(kk)
if(kk.eq.nstrt)goto 341
goto 361
endif
goto 351
endif
enddo
351 enddo
enddo
361 enddo


341 jj=0
do i=1,kk
jj=jj+1
do m18=1,nae
str3(jj,m18)=str1(st_ct1(i),m18)
enddo
q_fac2(jj)=q_fac1(st_ct1(i))
stq3(jj)=stq1(st_ct1(i))
stq4(jj)=stq2(st_ct1(i))
bdq1(jj)=bdq(st_ct1(i))
enddo

do i=1,n
do i2=1,kk
if(i.eq.st_ct1(i2))goto 371
enddo
jj=jj+1
do m18=1,nae
str3(jj,m18)=str1(i,m18)
enddo
q_fac2(jj)=q_fac1(i)
stq3(jj)=stq1(i)
stq4(jj)=stq2(i)
bdq1(jj)=bdq(i)
371 enddo

print*,'souravooooo'
do i=1,n
q_fac1(i)=0
str_quality_1(i)=0
str_quality_2(i)=0
bondq(i)=0
do i2=1,nae
str1(i,i2)=0
enddo
enddo

do i=1,n
q_fac1(i)=q_fac2(i)
str_quality_1(i)=stq3(i)
str_quality_2(i)=stq4(i)
bondq(i)=bdq1(i)
do i2=1,nae
str1(i,i2)=str3(i,i2)
enddo
enddo

endif


endif

do m18=1,n
Print*,'str1*******',m18,(str1(m18,m19),m19=1,nae),q_fac1(m18),bondq(m18)
!,str_quality_1(m18),str_quality_2(m18),bondq(m18),qulsym(m18),symq(m18)
!print*,'q_fac1',q_fac1(m18)
enddo
!do i=1,kk
!print*,'st_ct1(kk)',st_ct1(i),kk
!enddo
!print*,'i4',i4
deallocate(str3)
deallocate(q_fac3)
deallocate(q_fac5)
deallocate(q_fac6)
deallocate(q_fac7)
deallocate(q_fac2)
deallocate(qul1)
deallocate(qual3)
deallocate(qual4)
deallocate(bdq)
deallocate(bdq1)
deallocate(str_cnt)
deallocate(stq1)
deallocate(stq2)
deallocate(stq3)
deallocate(stq4)
deallocate(sl_num)
deallocate(sl_qfac)
deallocate(sl_qfac1)
print*,'exit qult_str_arrange'

!stop
return
end subroutine qult_str_arrange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine prio_rad_str(nl,str1,ncqs,pref_radical)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,j,i3,i4,i5,i6,i7,ii,m16,m17,m18,m19,m23,nl,ncqs,qual2(15000),qul1(15000),&
qul2(15000),jj,nqul1,nqul2,qt,qual3(15000),qual4(15000),k6,k7,k8,pref_radical(15000)
integer::str1(15000,20),qual1(15000),qcs_serial(15000),q_fac(15000)

!!!!! starting the arrangement of the structures according to priority radicles!!!


print*,'enter prio_rad_str'
do i=1,15000
!pref_radical(i)=2
pref_radical(i)=nlast+1
enddo

!print*,'nlpset',nlpset
if(nlpset.eq.0)goto 102
do i5=1,nlpset
ii=0
do i=1,nl*2,2
!print*,'i',i
do j=1,nl
if(plpair(i5,j).eq.str1(1,i))then
ii=ii+1
endif
enddo
enddo
!print*,'ii',ii
if(ii.ne.nl) goto 103
if(ii.eq.nl)then
jj=i5
if(plpair(i5,nl+1).ne.0.)goto 100
if(plpair(i5,nl+1).eq.0.)goto 102
endif
103 enddo
goto 101
!102 print*,'i51',jj,plpair(i5,nl+1)


!102 print*,'sourav'
102 do i=1,ncqs
!print*,'sourav1',(str1(i,i4),i4=1,nae)
do j=1,prad
!print*,'norad',norad(j)
ii=0
do i3=1,norad(j)
do i4=nae-nlast+1,nae
if(str1(i,i4).eq.prio_rad(j,i3))then
ii=ii+1
endif
enddo
enddo
if(ii.eq.norad(j))then
pref_radical(i)=pref_radical(i)-1
!print*,'pref_radical(i)',pref_radical(i)
endif
enddo
enddo

!print*,'sourav2'
if(nlpset.eq.0)goto 101
100 if(plpair(jj,nl+1).eq.0) goto 101


do i=1,ncqs
do j=nl+1,lp(jj)
i7=plpair(jj,j)
ii=0
do i3=1,norad(i7)
do i4=nae-nlast+1,nae
if(str1(i,i4).eq.prio_rad(i7,i3))then
ii=ii+1
endif
enddo
enddo
!print*,'norad(i7)',norad(i7)
if(ii.eq.norad(i7))then
pref_radical(i)=pref_radical(i)-1
endif
enddo
enddo

!101 enddo


!101 do m19=1,ncqs
!!print*,'pref_rad',pref_radical(m19)
!write(*,231),m19,(str1(m19,i4),i4=1,nae)
!enddo
231 format(30I3)

print*,'exit prio_rad_str'

101 return
end subroutine prio_rad_str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine main_bond_str(nl,str1,ncqs,qual1,qual2,str2,q_fac)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::jj,nl,str1(15000,15),str2(15000,15),ncqs,qual1(15000),qual2(15000),i,j,k,l &
,qul2(15000),qul1(15000),q_fac(15000),nqul,m19,k6,k7,k8,nqul1,nqul2,qt,qtg &
,qual3(15000),qual4(15000)

print*,'enter main_bond_str'
!do i=1,ncqs
!write(*,231)qual1(i),qual2(i),(str1(i,j),j=1,nae)
!enddo

jj=1
do m19=1,ncqs
!print*,'q_fac2',q_fac(m19)
if(m19.eq.1)qul1(1)=qual1(1)
j=jj
do i=1,j
if(qul1(i).eq.qual1(m19))goto 373
enddo
jj=jj+1
qul1(i)=qual1(m19)
!print*,'qul',qul(i)
373 enddo
nqul1=jj

jj=1
do m19=1,ncqs
!print*,'q_fac2',q_fac(m19)
if(m19.eq.1)qul2(1)=qual2(1)
j=jj
do i=1,j
if(qul2(i).eq.qual2(m19))goto 374
enddo
jj=jj+1
qul2(i)=qual2(m19)
!print*,'qul',qul(i)
374 enddo
nqul2=jj

do k6=1,15000
qual3(k6)=0
qual4(k6)=0
enddo

do k6=1,nqul1
do k7=1,nqul1
do k8=1,k6
if(qual3(k8).eq.qul1(k7))goto 389
enddo
if(qul1(k7).gt.qual3(k6))then
qual3(k6)=qul1(k7)
!print*,'**',qul1(k6)
endif
389 enddo
enddo

do k6=1,nqul2
do k7=1,nqul2
do k8=1,k6
if(qual4(k8).eq.qul2(k7))goto 399
enddo
if(qul2(k7).gt.qual4(k6))then
qual4(k6)=qul2(k7)
!print*,'**',qul1(k6)
endif
399 enddo
enddo

do k6=1,15000
qul1(k6)=0
qul2(k6)=0
enddo

k7=0
do k6=nqul1,1,-1
k7=k7+1
qul1(k7)=qual3(k6)
enddo

k7=0
do k6=nqul2,1,-1
k7=k7+1
qul2(k7)=qual4(k6)
enddo

k=0
qt=0
qtg=1
do m19=1,nqul1
!do j=((nae-(nl*2+nlast))/2)+1,1,-1
do j=1,nqul2
if(qtg.ne.k)qt=qt+1
qtg=k
!print*,'jjj',j
do i=1,ncqs
!print*,'qul,q_fac',qul(m19),q_fac(i),nqul
if(qul1(m19).ne.qual1(i))goto 100
!if(qual2(i).eq.j-1)then
if(qual2(i).eq.qul2(j))then
k=k+1
do l=1,nae
str2(k,l)=str1(i,l)
enddo
!write(*,231)(str2(k,k7),k7=1,nae)
q_fac(k)=qt
!print*,qt
endif
100 enddo
enddo
enddo

!do i=1,nqs
!write(*,231)(str2(i,j),j=1,nae)
!enddo
231 format(20I3)
print*,'exit main_bond_str'

return
end subroutine main_bond_str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine main_bond_cal(nl,str1,ncqs,mbondq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::m19,m18,i1,i2,i3,i4,i5,i6,i7,i8,i9,ii,iii,iiii,nl,ncqs
integer::nn(10),str1(15000,20),mbondq(15000)

print*,'enter main_bond_cal'
do i1=1,15000
mbondq(i1)=1+(nae-nl*2-nlast)/2
enddo

!!do i1=1,6
!!write(*,231)(main_bond(i1,i2),i2=1,nmbond(i1))
!!enddo
!i4=0
!do i1=1,ncqs
!!write(*,231)(str1(i1,i2),i2=1,nae)
!i8=1
!if(nl.ne.0)then
!do i3=1,nlpset
!ii=0
!do i6=1,nl*2
!do i4=1,lp
!!print*,'plpair',plpair(i3,i4)
!if(str1(i1,i6).eq.plpair(i3,i4))then
!ii=ii+1
!endif
!enddo
!enddo
!if(ii.eq.nl*2)then
!i8=i3
!!print*,'iiii88888',i8
!goto 529
!endif
!678 enddo
!goto 527
!endif


do i1=1,ncqs
i9=0
do i3=1,100
bond_count(i1,i3)=0
enddo
iiii=1+(nae-nl*2-nlast)/2
do i3=1+nl*2,(nae-nlast),2
i9=0
nn(1)=0
nn(2)=0
do i4=1,nmbond*2,2
i9=i9+1
iii=0
do i5=i3,i3+1
do i7=i4,i4+1
!print*,str1(i1,i5),main_bond(i8,i7),i7,i8
if(str1(i1,i5).eq.main_bond(i7))then
!print*,'***********',str1(i1,i5),main_bond(i8,i7)
iii=iii+1
!print*,'iii',iii
if(mod(i5,2).eq.0) nn(2)=i5
if(mod(i5,2).eq.1) nn(1)=i5
endif
enddo

enddo
if(iii.ne.2) goto 517
iiii=iiii-1
bond_count(i1,i9)=bond_count(i1,i9)+1
517 enddo
enddo
mbondq(i1)=iiii
!print*,'iiii',iiii
enddo

do i1=1,ncqs
write(*,231)(str1(i1,i2),i2=1,nae),mbondq(i1)
write(*,231)(bond_count(i1,i2),i2=1,nmbond)
enddo
231 format (30I3)

print*,'exit main_bond_cal'
return
end subroutine main_bond_cal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine nnat_bond_loop(k1,k2,new_mat,nrow)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!implicit none 
!
!integer::i,i1,j,k1,k1,nrow,new_mat(10,10)
!
!do i=1,k1
!do i1=1,k2
!
!enddo
!enddo
!
!end subroutine nnat_bond_loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symmetry_cal_pi(nl,str1,ncqs,stsymsc,symq,nssym)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

common/ats/atsymset,nsym,syn,at_sym
common/loops/full_nn_group,fullgrp,natom,nelimt,sl_group,tot_orb,nn_group

integer::m19,m18,i,i1,i2,i3,i4,i5,i6,i7,i8,ii,iii,iiii,nl,ncqs,k,k1,k2,n1,n2,n3,n4,n5,n6,l6,&
nd,fullgrp,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,&
k25,k26,k27,n,j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,nssym,j,jj,nn1,nn2,kkk,n11,n12,n17,n18,n9,l,b3,b4
integer::nnat_bond_new(100,2),tot_orb(100),nn_count(50),sig_orb(100),sig_orb_1(100),nsig,nnatom,natom
integer::kk,nscore(100),nnscore,ipnum,connectivity_row(20,20),b1,b2,l1,l2,l3,l4,nnat_bond_1(100,2)
double precision::tscore,stsymsc(15000),strscore,bd_dist,score(100,2),new_score(100,2),lpna,&
ssym(15000),order(15000)
integer::nn(10),str1(15000,20),mbondq(15000),nn_group(50,10),nelimt(50),full_nn_group(1000),&
sl_group(50,10),symq(15000),numbond,loop_score,nrow,new_mat(10,10),ncrow

integer::atsymset(20,20),nsym,syn(50),at_sym(50),loop_score_row(20),coordination_mat(20,20)
real*8::atoset_symscr(1000),piscr,losc,lone_score(15000),symq1(15000),cot(1000),sig_sym_sc(15000)&
,sym_sets(15000),n7,n8,coord_score,bndscore(20)
real::stsymsc1(15000)

print*,'symm_pi'
open(unit=15,file='sympi.temp',status='unknown')

do i3=1,atom
print*,'atoset_symm_pi',(atoset(active_atoms(i3),i4),i4=1,atn(active_atoms(i3)))
enddo


do i=1,100
  do i1=1,2
    score(i,i1)=0.0
  enddo
enddo

i4=0
do i3=1,atom
  lpna=0.0
  k=i3
  if(input_flg.eq.1)then
    do i1=1,niabd
      if(k.eq.iabd(i1))then
       lpna=lpna+1.0/prime_num(2)
      print*,'lpna1',lpna
      endif
    enddo
    do i2=1,nialp
      if(k.eq.ialp(i2))then
       lpna=lpna+1.0/prime_num(1)
      print*,'lpna2',lpna
      endif
    enddo
    do i2=1,niach
      if(k.eq.iach(i2))then
       lpna=lpna+1.0/prime_num(3)
      print*,'lpna3',lpna
      endif
    enddo
  else
!  lpna=lpna+biasval(i3)+dist_nnat(i3)+atn(active_atoms(i3))/prime_num(1)
  lpna=biasval(i3)+dist_nnat(i3)
      print*,'lpna4',lpna
  endif
!   score(i3,2)=lpna+1.0/at_num(i3)
   score(i3,2)=lpna+at_num(i3)
!print*,i3,lpna,at_num(i3),score(i3,2)
enddo

!!!!!!!!!!!!!!!! coordination score start !!!!!!!!!!!!!!
do i3=1,atom
print*,'atoset',(atoset(i3,i4),i4=1,atn(i3))
enddo
do i1=1,nnnatom
print*,'nnat_bond1',(nnat_bond(i1,i2),i2=1,2)
enddo

!nnatom=0
!do i1=1,atom
!do i3=1,atn(nnat_bond(i1,1))
!do i4=1,atn(nnat_bond(i1,2))
!nnatom=nnatom+1
!nnat_bond_new(nnatom,1)=atoset(nnat_bond(i1,1),i3)
!nnat_bond_new(nnatom,2)=atoset(nnat_bond(i1,2),i4)
!enddo
!enddo
!enddo

!do i1=1,nnatom
!print*,'nnat_bond_new',(nnat_bond_new(i1,i2),i2=1,2)
!enddo
!stop

!l1=100000
!l2=0
!do i1=1,nnnatom
!do i2=1,2
!if(nnat_bond(i1,i2).lt.l1)l1=nnat_bond(i1,i2)
!if(nnat_bond(i1,i2).gt.l2)l2=nnat_bond(i1,i2)
!enddo
!enddo
!print*,'l1,l2',l1,l2,l2-l1+1,nnnatom*2

!if(l2-l1+1.gt.nnnatom*2)then
!l3=0
!l4=0
!do i1=l1,l2
!if(l4.eq.0)l3=l3+1
!l4=1
!do i2=1,nnnatom
!do i3=1,2
!!print*,'i1,i2,i3,l3',i1,i2,i3,l3
!if(nnat_bond(i2,i3).eq.i1)then
!nnat_bond_1(i2,i3)=l3
!!print*,'i1,i2,i3,l3**',i1,i2,i3,l3
!l4=0
!endif
!enddo
!enddo
!enddo
!!endif
!k=0
!do i1=1,nnnatom
!do i2=1,nrow
!k=k+1
!do i3=1,2
!nnat_bond_1(k,i3)=new_mat(i2,i3)
!enddo
!enddo
!enddo
!
!do i1=1,nnnatom
!print*,'nnat_bond2',(nnat_bond_1(i1,i2),i2=1,2)
!enddo



!n6=0
!n7=0
!n8=0
!do i2=1,nnnatom
!do i1=1,2
!nnat_bond_new(i2,i1)=nnat_bond(i2,i1)
!enddo
!print*,(nnat_bond_new(i2,i1),i1=1,2),'|',(nnat_bond_1(i2,i1),i1=1,2)
!enddo

n5=0
!do i2=1,nnatom
do i2=1,nnnatom
do i1=1,2
do i3=1,n5
!if(nnat_bond_new(i2,i1).eq.tot_orb(i3)) goto 531
if(nnat_bond(i2,i1).eq.tot_orb(i3)) goto 531
enddo
n5=n5+1
!tot_orb(n5)=nnat_bond_new(i2,i1)
tot_orb(n5)=nnat_bond(i2,i1)
531 enddo
enddo



!do i1=1,atom
!n5=n5+1
!tot_orb(n5)=i1+nao+niao
!enddo

do i2=1,n5
print*,'n5,tot_orb',n5,tot_orb(i2)
enddo

do i1=1,atom
n1=0
!do i3=1,nnatom
do i3=1,nnnatom
do i4=1,2
!if(active_atoms(i1).eq.nnat_bond_new(i3,i4))then
if(active_atoms(i1).eq.nnat_bond(i3,i4))then
n1=n1+1
endif
enddo
enddo
nn_count(i1)=n1
enddo
do i3=1,atom
print*,'nn_count',nn_count(i3)
enddo

print*,'atom',atom,nnnatom

m18=0
m19=0
do i3=1,n5
k=tot_orb(i3)
m19=0
!do i2=1,nnatom
do i2=1,nnnatom
do i1=1,2
!print*,'atom',k,nnat_bond_new(i2,i1)
!if(k.eq.nnat_bond_new(i2,i1))then
if(k.eq.nnat_bond(i2,i1))then
if(i1.eq.1)i5=2
if(i1.eq.2)i5=1
m19=m19+1
m18=m18+1
!nn_group(i3,m19)=nnat_bond_new(i2,i5)
nn_group(i3,m19)=nnat_bond(i2,i5)
print*,'nnat_bond_new(i2,i5)',i2,i5,nnat_bond(i2,i5)
sl_group(i3,m19)=m18
endif
enddo
enddo
nelimt(i3)=m19
!print*,'nelimt(i3)',nelimt(i3)
enddo

!!! 'nn_group' specify orbitals in tot_orb array are connected with which other
!!! orbitals

natom=n5
do i3=1,n5
print*,'*nn_group*',nelimt(i3),(nn_group(i3,i4),i4=1,nelimt(i3))
enddo
do i3=1,n5
print*,'*sl_group*',(sl_group(i3,i4),i4=1,nelimt(i3))
enddo

m19=0
do i3=1,n5
do i4=1,nelimt(i3)
m19=m19+1
full_nn_group(m19)=nn_group(i3,i4)
enddo
enddo
fullgrp=m19

write(*,*),'full_nn_g',m19,(full_nn_group(i3),i3=1,m19)


!!!!!!!!!!!!!!!! coordination score ends !!!!!!!!!!!!!!
do i=1,ncqs
do i3=1,atom
new_score(i3,2)=score(i3,2)
print*,'score',i,i3,score(i3,2)
enddo
print*,'*****************************************************'
write(*,233),'structures',(str1(i,i3),i3=1,nae)
233 format(a,20I3)

!! if the structures have active lone pairs and or radicals the associated atoms are being scored in 'new_score()'.
!!if(nl.ne.0)then
!!  do i2=1,nl*2,2
!!    do i3=1,atom
!!      do i4=1,atn(active_atoms(i3))
!!        if(str1(i,i2).eq.atoset(active_atoms(i3),i4))then
!!          new_score(i3,2)=new_score(i3,2)+1.0/prime_num(2)
!!        endif
!!      enddo
!!    enddo
!!  enddo
!!endif
!!
!!if(nlast.ne.0)then
!!  do i2=nae-nlast+1,nae
!!    do i3=1,atom
!!      do i4=1,atn(active_atoms(i3))
!!        if(str1(i,i2).eq.atoset(active_atoms(i3),i4))then
!!          new_score(i3,2)=new_score(i3,2)+1.0/prime_num(3)
!!        endif
!!      enddo
!!    enddo
!!  enddo
!!endif


!! Scoring of the bonds of the structures started from here
!! the score of the bond is given by the score of the atoms involeved in the bond

strscore=0.0
!ipnum=3
ipnum=0
numbond=0
do i1=1+nl*2,nae-nlast,2
tscore=0.0
n1=0
n2=0
n3=0
n4=0
n5=0
n6=0
n11=0
n12=0
n17=0
n18=0
n9=0
kkk=0
  do i2=i1,i1+1
  b1=str1(i,i1)
  b2=str1(i,i1+1)
print*,'bond between',str1(i,i1),str1(i,i1+1)
    do i3=1,atom
      do i4=1,atn(active_atoms(i3))
        if(str1(i,i2).eq.atoset(active_atoms(i3),i4))then
          if(i2.eq.i1) then
          n1=i3
          n11=active_atoms(i3)
            do k1=1,nsym
              do k2=1,syn(k1)
                if(atoset(active_atoms(i3),i4).eq.atsymset(k1,k2))then
                n3=k1
                n18=atsymset(k1,k2)
                endif
              enddo
            enddo
          endif
          if(i2.eq.i1+1) then
          n2=i3
          n12=active_atoms(i3)
            do k1=1,nsym
              do k2=1,syn(k1)
                if(atoset(active_atoms(i3),i4).eq.atsymset(k1,k2))then
                n4=k1
                n5=atsymset(k1,k2)
                endif
              enddo
            enddo
          endif
        endif
      enddo
    enddo
  enddo

!print*,'n1,n2,n3,n4',n1,n2,n3,n4
if(n3.eq.n4)then
  do k1=1,nsym
    do k2=1,syn(k1)
      if(n5.eq.atsymset(k1,k2))then
      n6=at_sym(k1)
      endif
    enddo
  enddo
!print*,'if orb similar pi or sig',n6
if(n6.eq.1)kkk=ipnum+2
if(n6.eq.2)kkk=ipnum+1
!  tscore=new_score(n1,2)+new_score(n2,2)+1.0/prime_num(kkk)
  tscore=new_score(n1,2)+new_score(n2,2)+kkk/(new_score(n1,2)+new_score(n2,2))
!print*,'score_similar',tscore,new_score(n1,2),new_score(n2,2),kkk/(new_score(n1,2)+new_score(n2,2))
endif


if(n3.ne.n4)then
  do k1=1,nsym
    do k2=1,syn(k1)
      if(n5.eq.atsymset(k1,k2))then
      n6=at_sym(k1)
      endif
    enddo
  enddo
do k1=1,nsym
  do k2=1,syn(k1)
    if(n18.eq.atsymset(k1,k2))then
     n9=at_sym(k1)
    endif
  enddo
enddo

!print*,'if orb desimilar',n6,n9
if(symtype.eq.'loose')then
  if(n6.eq.1.and.n9.eq.1)kkk=ipnum+2
    if(n6.ne.n9)kkk=ipnum+3
  endif
!print*,'kkk',kkk
if(symtype.eq.'tight')then
  if(n6.eq.1.and.n9.eq.1)then
  kkk=ipnum+3
!print*,'kkk3',kkk
  else
  if(n3.eq.1.or.n4.eq.1)kkk=ipnum+4
  if(n3.eq.2.or.n4.eq.2)kkk=ipnum+5
!print*,'kkk4,5',kkk
  endif
endif
!tscore=new_score(n1,2)+new_score(n2,2)+1.0/prime_num(kkk)
tscore=new_score(n1,2)+new_score(n2,2)+kkk/(new_score(n1,2)+new_score(n2,2))
!print*,'score_disimilar',tscore,new_score(n1,2),new_score(n2,2),kkk
endif
nnscore=0.0
nd=0

do i2=1,nactorb
  if(str1(i,i1).eq.active_orbs(i2))then
  k1=atm_nb_orbs(i2)
  endif
  if(str1(i,i1+1).eq.active_orbs(i2))then
  k2=atm_nb_orbs(i2)
  endif
enddo
tscore=tscore+dist_mat(k1,k2)
!print*,'tscore with bond length',tscore,dist_mat(k1,k2)
print*,'b1,b2****',b1,b2,atom
do i3=1,atom
do i4=1,atn(active_atoms(i3))
print*,'i3,i4,b1,atoset(i3,i4)',i3,i4,b1,atoset(active_atoms(i3),i4)
if(b1.eq.atoset(active_atoms(i3),i4))then
b3=active_atoms(i3)
endif
enddo
enddo

do i3=1,atom
do i4=1,atn(active_atoms(i3))
if(b2.eq.atoset(active_atoms(i3),i4))then
b4=active_atoms(i3)
endif
enddo
enddo

print*,'b3,b4*****',b3,b4

do i2=1,natom
if(b3.eq.tot_orb(i2))l1=i2
enddo
!print*,'l1',l1
do i2=1,nelimt(l1)
 if(b4.eq.nn_group(l1,i2))then
!n=1
!nscore(n)=prime_num(9)
!nscore(n)=0
!print*,'nscore before loop start',nscore(n)
!goto 101
loop_score=1
numbond=numbond+1
loop_score_row(numbond)=loop_score+1
coordination_mat(numbond,1)=b3
coordination_mat(numbond,2)=b4
bndscore(numbond)=tscore
goto 555
 endif
enddo

 call symm_loops(b3,b4,loop_score,connectivity_row,ncrow)
print*,'ncrow',ncrow
do i3=1,ncrow
print*,'loop_score 1',loop_score
print*,'connectivity_row',(connectivity_row(i3,i2),i2=1,loop_score)
enddo
do i3=1,ncrow
numbond=numbond+1
loop_score_row(numbond)=loop_score+1
do i2=1,loop_score+1
coordination_mat(numbond,i2)=connectivity_row(i3,i2)
enddo
bndscore(numbond)=tscore
enddo

!bndscore(numbond)=tscore
555 strscore=strscore+tscore
print*,'strscore',strscore
enddo
call coordination_val(coordination_mat,numbond,loop_score_row,bndscore,coord_score)
print*,'str_score',strscore,coord_score
stsymsc(i)=1.0/strscore+coord_score
write(15,*)stsymsc(i)
enddo

rewind(15)
do i=1,ncqs
read(15,105)stsymsc1(i)
enddo
105 format(F8.4)

jj=1
do m19=1,ncqs
if(m19.eq.1)ssym(1)=stsymsc1(1)
j=jj
do i=1,j
if(ssym(i).eq.stsymsc1(m19))goto 373
enddo
jj=jj+1
ssym(i)=stsymsc1(m19)
373 enddo
nssym=jj
print*,'nssym',nssym

do i=1,10000
order(i)=0.0
enddo

do k3=1,nssym
do k4=1,nssym
do k5=1,k3
if(order(k5).eq.ssym(k4))goto 389
enddo
if(ssym(k4).gt.order(k3))then
order(k3)=ssym(k4)
endif
389 enddo
enddo

n=0
jj=0
do i=1,nssym
jj=jj+1
do j=1,ncqs
if(order(i).eq.stsymsc1(j))then
symq(j)=jj
endif
enddo
enddo
print*,'symqsymq_pi'
do j=1,ncqs
write(*,231),j,(str1(j,k),k=1,nae),symq(j)
enddo

231 format(20I3)

CALL SYSTEM ("rm sympi.temp")
200 return
end subroutine symmetry_cal_pi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symmetry_cal_sig(nl,str1,ncqs,stsymsc,symq,nssym)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

integer::m19,m18,i,i1,i2,i3,i4,i5,i6,i7,i8,ii,iii,iiii,nl,ncqs,k,k1,k2,n1,n2,n3,n4,n5,n6,l1,l3,l6,&
nd,fullgrp,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,&
k25,k26,k27,n,j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,nssym,j,jj,nn1,nn2,kkk,n17,n18,n9,ii2,nnat_bond_new_1(100,2)
integer::nnat_bond_new(100,2),tot_orb(100),nn_count(50),sig_orb(100),sig_orb_1(100),nsig,nnatom,natom&
,nn4,nn5,nna
real*8::tscore,score(100,2),new_score(100,2),stsymsc(15000),ssym(15000),order(15000),lpna,&
losc,lone_score(15000),symq1(5000),cot(1000),sig_sym_sc(15000),sym_sets(15000),n7,n8,bd_dist
integer::kk,nscore(100),nnscore,ipnum,loop_score,loopsc,ncrow
!double precision::rscore,strscore
real*8::rscore,strscore,scr,bndscore(20),coord_score
!real::nnscore,rscore,strscore
integer::nn(10),str1(15000,20),mbondq(15000),nn_group(50,10),nelimt(50),&
full_nn_group(1000),sl_group(50,10),symq(15000),connectivity_row(20,20)
!integer::kval(1000),deadend(20)
integer::atsymset(20,20),nsym,syn(50),at_sym(50),numbond,numbond1,loop_score_row(20),coordination_mat(20,20)
real*8::atoset_symscr(1000),piscr
real::stsymsc1(15000)

common/ats/atsymset,nsym,syn,at_sym
common/loops/full_nn_group,fullgrp,natom,nelimt,sl_group,tot_orb,nn_group
print*,'enter symmetry_cal_sig'
!! iabd--> inactive bonds associated with the atoms
!! ialp--> inactive lone pairs associated with the atoms
!! score()--> scoring of the atoms depending on the the 'iabd', 'ialp' and atomic number (at_number)
open(unit=15,file='symsig.temp',status='unknown')

do k1=1,nsym
print*,'at_sym(k1)',at_sym(k1)
enddo



!do i3=1,atom
!print*,'atoset_symm_sig',(atoset(active_atoms(i3),i4),i4=1,atn(active_atoms(i3)))
!print*,'active_atoms',active_atoms(i3)
!enddo


do i=1,100
do i1=1,2
score(i,i1)=0.0
enddo
enddo

print*,'prime_num',prime_num(2)
i4=0
do i3=1,atom
!lpna=prime_num(1)
lpna=0.0
k=i3
!k=active_atom(i3)
if(input_flg.eq.1)then
do i1=1,niabd
if(k.eq.iabd(i1))then
lpna=lpna+1.0/prime_num(2)
print*,'lpna',lpna,prime_num(2)
endif
enddo
!print*,'rajat roy',nialp
do i2=1,nialp
if(k.eq.ialp(i2))then
!print*,'rajat'
lpna=lpna+1/prime_num(1)
endif
enddo
do i2=1,niach
if(k.eq.iach(i2))then
!print*,'rajat'
lpna=lpna+1/prime_num(3)
endif
enddo
else

!lpna=lpna+biasval(i3)+dist_nnat(i3)+atn(active_atoms(i3))/prime_num(1)
lpna=lpna+biasval(i3)+dist_nnat(i3)
endif
!score(i3,2)=lpna+1.0/prime_num(at_num(i3))
!score(i3,2)=lpna+1.0/at_num(i3)
score(i3,2)=lpna+at_num(i3)
print*,'score:symme',score(i3,2),at_num(i3),biasval(i3),dist_nnat(i3)
enddo



print*,'symmetry_cal',nnnatom
n6=0
n7=0
n8=0
do i2=1,nnnatom
do i1=1,2
nnat_bond_new(i2,i1)=nnat_bond(i2,i1)+nao+niao
enddo
print*,(nnat_bond_new(i2,i1),i1=1,2),'|',(nnat_bond(i2,i1),i1=1,2)
enddo
n5=0
do i2=1,nnnatom
do i1=1,2
do i3=1,n5
if(nnat_bond_new(i2,i1).eq.tot_orb(i3)) goto 537
enddo
n5=n5+1
tot_orb(n5)=nnat_bond_new(i2,i1)
537 enddo
enddo
!n5=0
!do i1=1,atom
!n5=n5+1
!tot_orb(n5)=i1+nao+niao
!enddo

do i1=1,atom
n1=0
do i3=1,nnnatom
do i4=1,2
if(active_atoms(i1).eq.nnat_bond(i3,i4))then
!if(i1.eq.nnat_bond(i3,i4))then
n1=n1+1
endif
enddo
enddo
nn_count(i1)=n1
enddo
do i3=1,atom
print*,'nn_count',nn_count(i3)
enddo



! do i3=1,atom
!print*,'active_atom',i3,active_atoms(i3)
! do i4=1,atn(active_atoms(i3))
!print*,'atn',i4,atn(active_atoms(i3))
!print*,'atn(i3),syn(2)',syn(nsym),(atsymset(nsym,k2),k2=1,syn(nsym))
!enddo
!enddo

do i3=1,atom
print*,'atoset',(atoset(i3,i4),i4=1,atn(i3))
enddo


print*,'atom:atom',atom
n2=nnnatom
n4=0
n3=0
n6=0
 do i3=1,atom
print*,'iiiiiiiiiiiiiiiiiiii33333333333',i3
n1=0
do i4=1,atn(active_atoms(i3))
! do i4=1,atn(i3)
!print*,'atn**',i4,atn(active_atoms(i3))
 do k2=1,syn(nsym)
!print*,'atn(i3),syn(2)',atn(active_atoms(i3)),syn(nsym),atsymset(nsym,k2)
 if(atoset(active_atoms(i3),i4).eq.atsymset(nsym,k2))then
! if(atoset(i3,i4).eq.atsymset(nsym,k2))then
n1=n1+1
n6=n6+1
!sig_orb(n1)=atoset(i3,i4)
sig_orb(n1)=atoset(active_atoms(i3),i4)
sig_orb_1(n6)=sig_orb(n1)
print*,'sig_orb(n1)',sig_orb(n1),n1
 endif
 enddo
 enddo
!!nsig=n1
 if(n1.gt.0)then
n3=n3+1
k=0
!if(nn_count(i3).eq.1)k=1
!print*,'k',k,i3,nn_count(i3)
 do i2=1,nnnatom
 do i1=1,2
 if(i3.eq.nnat_bond(i2,i1))then
k=k+1
if(k.gt.n1)goto 113
n4=n4+1
!sig_orb_1(n4)=sig_orb(k)
nnat_bond_new(i2,i1)=sig_orb(k)
print*,'nnat_bond_new(i2,i1)***',i2,i1,nnat_bond(i2,i1),nnat_bond_new(i2,i1),sig_orb(k)
 endif
 enddo
113 enddo
!print*,'n1',n1,niao+nao+n3,n3
 do i1=1,n1
n2=n2+1
n5=n5+1
nnat_bond_new(n2,1)=niao+nao+i3
nnat_bond_new(n2,2)=sig_orb(i1)
tot_orb(n5)=sig_orb(i1)
 enddo
if(n1.gt.2)then
 do i1=1,n1
if(i1.ne.n1)ii2=i1+1
if(i1.eq.n1)ii2=1
n2=n2+1
nnat_bond_new(n2,1)=sig_orb(i1)
nnat_bond_new(n2,2)=sig_orb(ii2)
 enddo
endif
endif
 enddo
!!!!!!!!!!!!
nnatom=n2
nsig=n6
!!!! nsig=number of sigma orbital
!!!! nnatom= number of connected pair of orbitals in the extended CM
print*,'nsig',nsig,nnatom

natom=n5
!!!! natom=number of extended orbitals

do i2=1,nnnatom
print*,'nnat_bond',(nnat_bond(i2,i3),i3=1,2)
enddo
do i2=1,nnatom
print*,'nnat_bond_new',(nnat_bond_new(i2,i3),i3=1,2),nnatom
enddo
!print*,'sig_orb_1(n4)',(sig_orb_1(i2),i2=1,n4)
do i2=1,n5
print*,'tot_orb',tot_orb(i2),n5,nnatom
enddo

call nnat_bond_sig(nnat_bond_new,nnatom,sig_orb_1,nsig,nnat_bond_new_1,nna)
do i2=1,nna
print*,'nnat_bond_new_new',(nnat_bond_new_1(i2,i3),i3=1,2),nna
enddo

!231 format(20I3)
m18=0
m19=0
do i3=1,n5
k=tot_orb(i3)
m19=0
!do i2=1,nnatom
do i2=1,nna
do i1=1,2
!print*,'atom',k,nnat_bond_new(i2,i1)
!if(k.eq.nnat_bond_new(i2,i1))then
if(k.eq.nnat_bond_new_1(i2,i1))then
if(i1.eq.1)i5=2
if(i1.eq.2)i5=1
m19=m19+1
m18=m18+1
nn_group(i3,m19)=nnat_bond_new_1(i2,i5)
!print*,'nnat_bond(i2,i5)',i3,k,i2,i1,i5,m19,nnat_bond_new(i2,i5),nn_group(i3,m19)
sl_group(i3,m19)=m18
endif
enddo
enddo
nelimt(i3)=m19
!print*,'nelimt(i3)',nelimt(i3)
enddo


!!! 'nn_group' specify orbitals in tot_orb array are connected with which other
!!! orbitals

do i3=1,n5
print*,'*nn_group*',nelimt(i3),(nn_group(i3,i4),i4=1,nelimt(i3))
enddo
do i3=1,n5
print*,'*sl_group*',(sl_group(i3,i4),i4=1,nelimt(i3))
enddo

m19=0
do i3=1,n5
do i4=1,nelimt(i3)
m19=m19+1
full_nn_group(m19)=nn_group(i3,i4)
enddo
enddo
fullgrp=m19

write(*,*),'full_nn_g',m19,(full_nn_group(i3),i3=1,m19)
!goto 200
!endif
!!!!! Scoring of the structures has been started from here !!!

do i=1,ncqs
sig_sym_sc(i)=1.0
enddo


do i=1,ncqs
loopsc=0
do i3=1,atom
new_score(i3,2)=score(i3,2)
enddo

!! if the structures have active lone pairs and or radicals the associated atoms are being scored in 'new_score()'.

!print*,'symmetry_cal_sym_cross'
write(*,231)(str1(i,i1),i1=1,nae)
!print*,'tscore**,nl',nl
!!if(nl.ne.0)then
!!!print*,'rajat',nl
!!do i2=1,nl*2,2
!!!print*,i2
!!do i3=1,atom
!!do i4=1,atn(active_atoms(i3))
!!!print*,'atoset(active_atoms(i3),i4)',atoset(active_atoms(i3),i4)
!!if(str1(i,i2).eq.atoset(active_atoms(i3),i4))then
!!
!!new_score(i3,2)=new_score(i3,2)+1.0/prime_num(2)
!!!print*,'new score2',new_score(i3,2)
!!endif
!!enddo
!!enddo
!!enddo
!!endif
!!
!!if(nlast.ne.0)then
!!!print*,'rajat',nlast,nae-nlast+1,atom
!!do i2=nae-nlast+1,nae
!!do i3=1,atom
!!do i4=1,atn(active_atoms(i3))
!!if(str1(i,i2).eq.atoset(active_atoms(i3),i4))then
!!
!!new_score(i3,2)=new_score(i3,2)+1.0/prime_num(3)
!!endif
!!enddo
!!!print*,'new score3',i3,new_score(i3,2)
!!enddo
!!enddo
!!endif


!! Scoring of the bonds of the structures started from here
!! the score of the bond is given by the score of the atoms involeved in the


!bond
!do i3=1,atom
!do i4=1,atn(active_atoms(i3))
!!print*,'atoset(active_atoms(i3),i4)',atoset(active_atoms(i3),i4)
!enddo
!enddo
! do i3=1,nsig
!  print*,'sig_orb_1(i3)',sig_orb_1(i3),nsig
!enddo

strscore=0.0
ipnum=3
!print*,'strset'
numbond=0
numbond1=0
do i1=1+nl*2,nae-nlast,2
tscore=0.0
rscore=0.0
n1=0
n2=0
n3=0
n4=0
do i2=i1,i1+1
do i3=1,atom
do i4=1,atn(active_atoms(i3))
!do i4=1,atn(i3)
if(str1(i,i2).eq.atoset(active_atoms(i3),i4))then
!if(str1(i,i2).eq.atoset(i3,i4))then
if(i2.eq.i1) then
n1=active_atoms(i3)
!n1=i3
endif
if(i2.eq.i1+1) then
n2=active_atoms(i3)
!n2=i3

endif
endif
enddo
enddo
!print*,'bond***',str1(i,i2),n1,n2
enddo
n4=str1(i,i1)
n5=str1(i,i1+1)
do i2=1,nsym
do i3=1,syn(i2)
!print*,'atsymset',i2,i3,n4,n4,atsymset(i2,i3)
if(n4.eq.atsymset(i2,i3))nn4=i2
if(n5.eq.atsymset(i2,i3))nn5=i2
enddo
enddo
!print*,'nn4,nn5',nn4,nn5,n4,n5
if(nsym.eq.1)then
if(symtype.eq.'tight'.or.symtype.eq.'loose')then
if(nn4.eq.1.and.nn5.eq.1)scr=1.0
endif
endif
if(nsym.eq.2)then
if(symtype.eq.'tight'.or.symtype.eq.'loose')then
if(nn4.eq.2.and.nn5.eq.2)scr=1.0
if(nn4.eq.1.and.nn5.eq.1)scr=2.0
if(nn4.eq.1.and.nn5.eq.2)scr=4.0
if(nn4.eq.2.and.nn5.eq.1)scr=4.0
endif
endif
if(nsym.eq.3)then
if(symtype.eq.'tight')then
if(nn4.eq.3.and.nn5.eq.3)scr=1.0
if(nn4.eq.1.and.nn5.eq.1)scr=2.0
if(nn4.eq.2.and.nn5.eq.2)scr=2.0
if(nn4.eq.1.and.nn5.eq.2)scr=3.0
if(nn4.eq.2.and.nn5.eq.1)scr=3.0
if(nn4.eq.1.and.nn5.eq.3)scr=4.0
if(nn4.eq.3.and.nn5.eq.1)scr=4.0
if(nn4.eq.2.and.nn5.eq.3)scr=5.0
if(nn4.eq.3.and.nn5.eq.2)scr=5.0
endif
if(symtype.eq.'loose')then
if(nn4.eq.3.and.nn5.eq.3)scr=1.0
if(nn4.eq.1.and.nn5.eq.1)scr=2.0
if(nn4.eq.2.and.nn5.eq.2)scr=2.0
if(nn4.eq.1.and.nn5.eq.2)scr=2.0
if(nn4.eq.2.and.nn5.eq.1)scr=2.0
if(nn4.eq.1.and.nn5.eq.3)scr=3.0
if(nn4.eq.3.and.nn5.eq.1)scr=3.0
if(nn4.eq.2.and.nn5.eq.3)scr=3.0
if(nn4.eq.3.and.nn5.eq.2)scr=3.0
endif
endif
!tscore=n4*new_score(n1,2)+n5*new_score(n2,2)
tscore=new_score(n1,2)+new_score(n2,2)+scr/(new_score(n1,2)+new_score(n2,2))
!print*,'tscore,scr,OSS,ABS(topology)',tscore,scr,scr/(new_score(n1,2)+new_score(n2,2)),new_score(n1,2)+new_score(n2,2)
nnscore=0.0
nd=0
n1=0
n2=0
n3=0
n4=0
!do i3=1,nsig
!print*,'sig_orb_1',sig_orb_1(i3)
!enddo
 do i2=i1,i1+1
!print*,'i2',i2,str1(i,i2)
 do i3=1,nsig
  if(sig_orb_1(i3).eq.str1(i,i2))then
n3=n3+1
if(i2.eq.i1)n1=str1(i,i2)
if(i2.eq.i1+1)n2=str1(i,i2)
!print*,'n1,n2',n1,n2,n3,str1(i,i2)
goto 232
  endif
 enddo
 do i3=1,atom
 do i4=1,atn(active_atoms(i3))
! do i4=1,atn(i3)
  if(str1(i,i2).eq.atoset(active_atoms(i3),i4))then
!  if(str1(i,i2).eq.atoset(i3,i4))then
   if(i2.eq.i1) then
n1=i3+nao+niao
   endif
    if(i2.eq.i1+1) then
n2=i3+nao+niao
    endif
  endif
 enddo
 enddo
 232 enddo

!print*,'n3,bons',n3,n1,n2,n5,n18
!print*,'tscore',tscore

!loop_score=1
if(n1.gt.n2)then
k1=n2
k2=n1
else
k1=n1
k2=n2
endif
print*,'k1,k2,bond',k1,k2,natom
do i2=1,natom
if(k1.eq.tot_orb(i2))l1=i2
enddo
!print*,'l1',l1
do i2=1,nelimt(l1)
 if(k2.eq.nn_group(l1,i2))then
n=1
!nscore(n)=prime_num(9)
!nscore(n)=0
!print*,'nscore before loop start',nscore(n)
!goto 101
ncrow=1
loop_score=1
numbond=numbond+1
loop_score_row(numbond)=loop_score+1
coordination_mat(numbond,1)=k1
coordination_mat(numbond,2)=k2
goto 555
 endif
enddo

if(k1.ne.k2)call symm_loops(k1,k2,loop_score,connectivity_row,ncrow)
do i3=1,ncrow
numbond=numbond+1
loop_score_row(numbond)=loop_score+1
do i2=1,loop_score+1
coordination_mat(numbond,i2)=connectivity_row(i3,i2)
enddo
print*,'coordination_mat',(coordination_mat(numbond,i2),i2=1,loop_score+1)
enddo
if(k1.eq.k2)loop_score=0
goto 555
!555 print*,'tscore_before',tscore
555 if(loop_score.ne.0)tscore=tscore+1.0/loop_score
!print*,'loop_score',loop_score
!print*,'kk*,tscore',tscore,loop_score
!endif
rscore=1.0/tscore
do i2=1,ncrow
numbond1=numbond1+1
bndscore(numbond1)=rscore
enddo
!print*,'nearest neighbor is',kk,tscore,rscore
strscore=strscore+rscore
!strscore=strscore+tscore
if(loop_score.ne.0)loopsc=loopsc+prime_num(loop_score)
if(loop_score.eq.0)loopsc=loopsc+1000
!print*,strscore
!print*,'strscore,tscore,loop_score,loopsc',strscore,tscore,prime_num(loop_score),loopsc
enddo
print*,'numbond',numbond
do i1=1,numbond
print*,'loop_score_row',loop_score_row(i1)
print*,'coordination_mat',(coordination_mat(i1,i2),i2=1,loop_score_row(i1))
enddo

call coordination_val(coordination_mat,numbond,loop_score_row,bndscore,coord_score)

stsymsc(i)=strscore+coord_score
!stsymsc(i)=strscore
loopsymsc(i)=loopsc
print*,'stsymsc(i)',i,stsymsc(i),loopsymsc(i),coord_score
write(15,*)stsymsc(i)
enddo



rewind(15)
do i=1,ncqs
read(15,909)stsymsc1(i)
print*,i,stsymsc1(i)
enddo
909 format (F10.6)


jj=1
do m19=1,ncqs
!print*,'q_fac2',q_fac2(m19)
if(m19.eq.1)ssym(1)=stsymsc1(1)
j=jj
do i=1,j
if(ssym(i).eq.stsymsc1(m19))goto 373
enddo
jj=jj+1
ssym(i)=stsymsc1(m19)
!print*,qul(i)
373 enddo
nssym=jj
print*,'nssym',nssym

do i=1,jj
print*,ssym(i)
enddo

do i=1,10000
order(i)=0.0
enddo

do k3=1,nssym
do k4=1,nssym
do k5=1,k3
if(order(k5).eq.ssym(k4))goto 389
enddo
if(ssym(k4).gt.order(k3))then
order(k3)=ssym(k4)
endif
389 enddo
print*,'**order',k3,order(k3)
enddo

jj=0
do i=1,nssym
jj=jj+1
do j=1,ncqs
if(order(i).eq.stsymsc1(j))then
write(*,231)(str1(j,i1),i1=1,nae)
print*,'order(i),stsymsc(j)',order(i),stsymsc1(j)
if(nlast.ne.0)then
losc=0.0
do i1=nae-nlast,nae
do i2=1,syn(1)
if(str1(j,i1).eq.atsymset(1,i2))goto 531
enddo
goto 532
531 losc=losc+1.0/prime_num(str1(j,i1))
532 enddo
lone_score(j)=losc
endif
symq(j)=jj
!lsym(j)=loopsymsc(j)
print*,'symq(j),lsym(j)',symq(j),loopsymsc(j)
endif
enddo
enddo




print*,'sub:symmetry_cal'
do i=1,ncqs
!if(symq(i).eq.4)then
write(*,231),(str1(i,j),j=1,nae),symq(i)
!print*,lone_score(i)
!endif
enddo
print*,'sub:symmetry_cal'
!do i=1,ncqs
!write(*,231)(str2(i,j),j=1,nae)
!enddo

231 format(20I3)
!call sigma_symm(nl,ncqs,str1)
CALL SYSTEM ("rm symsig.temp")

print*,'exit symmetry_cal_sig'
200 return
end subroutine symmetry_cal_sig
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nnat_bond_sig(nnat_bond_new,n,sig_orb,nsig,nnat_bond_new_1,j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,i1,i2,i3,i4,i5,i6,j,k,l,nnat_bond_new(100,2),nnat_bond_new_1(100,2),n
integer::atoset_sig(20,10),atn_sig(20),orbs_sig(10),sig_orb(100),nsig,k1,k2,c1,c2,s_orbs(0)
integer::k3,k4
real*8::coordx(100),coordy(100),coordz(100),a,b,c,d,sig_coord(50,3),norm,kk,sig_dist_mat(40,40)

common/coordinate/coordx,coordy,coordz


k=0
do i=1,atom
if(atn(i).gt.2)then
j=0
do i1=1,atn(i)
do i2=1,nsig
print*,'atoset',atoset(active_atoms(i),i1),sig_orb(i2)
if(atoset(active_atoms(i),i1).eq.sig_orb(i2))then
j=j+1
orbs_sig(j)=sig_orb(i2)
endif
enddo
enddo
if(j.gt.2)then
k=k+1
do i4=1,j
atoset_sig(k,i4)=orbs_sig(i4)
enddo
atn_sig(k)=j
endif
endif
enddo


do i=1,k
print*,'atoset_sig',(atoset_sig(i,j),j=1,atn_sig(i))
enddo
j=0
do i2=1,n
do i3=1,2
do i=1,k
do i1=1,atn_sig(i)
if(atoset_sig(i,i1).eq.nnat_bond_new(i2,i3))then
k1=nnat_bond_new(i2,i3)
if(i3.eq.1)k2=nnat_bond_new(i2,2)
if(i3.eq.2)k2=nnat_bond_new(i2,1)
do i4=1,atn_sig(i)
if(k2.eq.atoset_sig(i,i4))goto 103
enddo
endif
enddo
enddo
enddo
j=j+1
nnat_bond_new_1(j,1)=nnat_bond_new(i2,1)
nnat_bond_new_1(j,2)=nnat_bond_new(i2,2)
103 enddo



l=0
k1=0
do i=1,k
do i1=1,atn_sig(i)
do i2=1,n
do i3=1,2
if(atoset_sig(i,i1).eq.nnat_bond_new(i2,i3))then
k1=nnat_bond_new(i2,i3)
if(i3.eq.1)k2=nnat_bond_new(i2,2)
if(i3.eq.2)k2=nnat_bond_new(i2,1)
!endif
!enddo
if(k2.gt.nao+niao)k2=k2-nao-niao
do i4=1,atn_sig(i)
if(k2.eq.atoset_sig(i,i4))goto 100
enddo

c1=0
c2=0
do i5=1,atom
do i6=1,atn(i5)
if(k1.eq.atoset(active_atoms(i5),i6))c1=i5
if(k2.eq.atoset(active_atoms(i5),i6))c2=i5
enddo
enddo

if(c2.eq.0)c2=k2
if(c1.eq.c2)goto 100
l=l+1
print*,'c1, c2 ***************',c1,c2
a=coordx(c2)-coordx(c1)
b=coordy(c2)-coordy(c1)
c=coordz(c2)-coordz(c1)
norm=1.0/sqrt(a**2.0 + b**2.0 + c**2.0)
d=0.3
s_orbs(l)=k1
sig_coord(l,1)=coordx(c1)+a*d*norm
sig_coord(l,2)=coordy(c1)+b*d*norm
sig_coord(l,3)=coordz(c1)+c*d*norm

print*,'a,b,c',a,b,c
print*,'c1 coord',coordx(c1),coordy(c1),coordz(c1)
print*,'c2 coord',coordx(c2),coordy(c2),coordz(c2)
print*,'sig_coord',l,sig_coord(l,1),sig_coord(l,2),sig_coord(l,3)


endif
enddo
100 enddo
enddo

do i1=1,l
!print*,'s_orbs',s_orbs(i1)
print*,'sig_coord',i1,sig_coord(i1,1),sig_coord(i1,2),sig_coord(i1,3)
enddo
print*,'**************',sig_coord(2,1)

do i1=1,l
do i2=1,l
sig_dist_mat(i1,i2)=0.0
enddo
enddo

do i1=1,l
print*,'sig_coord i1** i1',i1,sig_coord(i1,1),sig_coord(i1,2),sig_coord(i1,3)
do i2=1,l
sig_dist_mat(i1,i2)=0.0
print*,'sig_coord',i2,sig_coord(i2,1),sig_coord(i2,2),sig_coord(i2,3)
sig_dist_mat(i1,i2)=sqrt((sig_coord(i1,1)-sig_coord(i2,1))**2.0 + &
(sig_coord(i1,2)-sig_coord(i2,2))**2.0 + &
(sig_coord(i1,3)-sig_coord(i2,3))**2.0)
!print*,'coord1',i1,sig_coord(i1,1),sig_coord(i1,2),sig_coord(i1,3)
!print*,'coord2',i2,sig_coord(i2,1),sig_coord(i2,2),sig_coord(i2,3)
print*,'sig_dist_mat',sig_dist_mat(i1,i2)
enddo
enddo

do i1=1,l
print*,'sig_dist_mat',(sig_dist_mat(i1,i2),i2=1,l)
enddo

do i1=1,l
kk=1000.0
do i2=1,l
if(i1.eq.i2)goto 101
if(sig_dist_mat(i1,i2).lt.kk)kk=sig_dist_mat(i1,i2)
101 enddo
print*,'kk',kk

do i2=1,l
if(kk.eq.sig_dist_mat(i1,i2))then


do i3=1,j
do i4=1,2
if(s_orbs(i1).eq.nnat_bond_new_1(i3,i4))then
k3=nnat_bond_new_1(i3,i4)
if(i4.eq.1)k4=nnat_bond_new_1(i3,2)
if(i4.eq.2)k4=nnat_bond_new_1(i3,1)
endif
enddo
if(s_orbs(i2).eq.k4) goto 102
enddo


j=j+1
nnat_bond_new_1(j,1)=s_orbs(i1)
nnat_bond_new_1(j,2)=s_orbs(i2)
print*,'nnat_bond_new_1',nnat_bond_new_1(j,1),nnat_bond_new_1(j,2)
endif
102 enddo
enddo

enddo

do i1=1,j
print*,'new_mat',(nnat_bond_new_1(i1,i2),i2=1,2)
enddo


return
end subroutine nnat_bond_sig
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine coordination_val(coordination_mat,numbond,loop_score_row,bndscore,coord_score)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,i1,i2,i3
integer::coordination_mat(20,20),numbond,loop_score_row(20),new_row,loop_score_new(20)
real*8::bndscore(20),coord_score,c_score,score
common/lp/loop_score_new

print*,'enter coordination_val ********************************'
do i=1,numbond
print*,'coordination_mat',(coordination_mat(i,i2),i2=1,loop_score_row(i))
enddo

do i=1,20
loop_score_new(i)=loop_score_row(i)
enddo

coord_score=0.0
do i=1,numbond
!print*,'coordination_mat,i',i,'>',(coordination_mat(i,i2),i2=1,loop_score_row(i))
if (coordination_mat(i,1).ne.0)then
c_score=0.0
do i1=i+1,numbond
!print*,'coordination_mat,i1',i1,'>',(coordination_mat(i1,i2),i2=1,loop_score_row(i1))
if(i.ne.i1)then
if (coordination_mat(i1,1).ne.0)then
score=0.0
do i2=1,loop_score_row(i)
do i3=1,loop_score_row(i1)
print*,'i2,i3',i2,i3
if (coordination_mat(i,i2).eq.coordination_mat(i1,i3))then
!score=score+(bndscore(i)/(loop_score_row(i)-1))*((loop_score_row(i)-1)-(new_row(i2,i)-1))+ &
!(bndscore(i1)/(loop_score_row(i1)-1))*((loop_score_row(i1)-1)-(new_row(i3,i1)-1))
score=score+(bndscore(i)/new_row(i2,i)+bndscore(i1)/new_row(i3,i1))
print*,'bndscore(i),bndscore(i1),(i,i2),(i1,i3),new_row(i2,i),new_row(i3,i1)',&
bndscore(i),bndscore(i1),'(',i,i2,')','(',i1,i3,')',new_row(i2,i),new_row(i3,i1)
!print*,'i,i1,i2,i3,score',i,i1,i2,i3,score
endif
enddo
enddo
if(score.ne.0.0)c_score=c_score+1.0/score
!c_score=c_score+score
!print*,'c_score',c_score
endif
endif
enddo
!if(c_score.ne.0.0)coord_score=coord_score+1.0/c_score
if(c_score.ne.0.0)coord_score=coord_score+c_score
print*,'coord_score',coord_score
print*,'*************'
endif
enddo


return
end subroutine coordination_val
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function new_row(l,n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer::n,l,i,i1,i2,i3,loop_score_row(20),new_row_1(20),new_row
common/lp/loop_score_row

k=int(loop_score_row(n)/2)+mod(loop_score_row(n),2)
if (mod(loop_score_row(n),2).eq.0)then
i1=0
do i=1,loop_score_row(n)/2
i1=i1+1
new_row_1(i1)=k-ABS(k-i)
enddo

do i=loop_score_row(n)/2,1,-1
i1=i1+1
new_row_1(i1)=k-ABS(k-i)
enddo

endif
if (mod(loop_score_row(n),2).eq.1)then
do i=1,loop_score_row(n)
new_row_1(i)=k-abs(k-i)
enddo
endif
new_row=new_row_1(l)
!print*,'new_row,n,new_row_1',new_row,n,'|',(new_row_1(i),i=1,loop_score_row(n))

return
end function new_row
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symm_loops(k1,k2,loop_score,connectivity_row,nrow)

use commondat
implicit none

common/loops/full_nn_group,fullgrp,natom,nelimt,sl_group,tot_orb,nn_group
integer::k1,k2,k3,k4,k5,k6,k7,i,i1,i2,i3,i4,i5,i6,j,natom,fullgrp,loop_score1,nrow
integer::nn_group(50,10),nelimt(50),full_nn_group(1000),sl_group(50,10),tot_orb(100)
integer::nloop(100),orbs1(200),orbs2(200),orbs3(40),gr_sl,n,nsl,nsl1,ii,iii,flg,loop_score
integer::ii0,ii1,ii2,ii3,ii4,ii5,ii6,ii7,ii8,ii9,ii10,ii11,ii12,ii13,j1,j2,&
j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,l,maxlp
integer::ls1,ls2,ls3,ls4,ls5,ls6,ls7,ls8,ls9,ls10,ls11,ls12,ls13,ls14,connectivity(20,20),&
connectivity_row(20,20),nn,nnn(10)

print*,'enter_symm_loops',k1,k2
!do i3=1,natom
!print*,'*nn_group*',nelimt(i3),(nn_group(i3,i4),i4=1,nelimt(i3))
!enddo
!do i3=1,natom
!print*,'*sl_group*',(sl_group(i3,i4),i4=1,nelimt(i3))
!enddo
nsl=0
nsl1=0
ii=0
iii=0
i4=0
i5=0
flg=0
maxlp=15
loop_score1=0
!do i=1,natom
!print*,'sourav1'
!print*,'tot_orb',tot_orb(i)
!enddo
!do i=1,natom
!print*,'nelimt',nelimt(i),tot_orb(i)
!enddo
l=0
do i=1,natom
if(k1.eq.tot_orb(i))goto 101
enddo

101 do j=1,nelimt(i)
!print*,'loop0'
do i6=1,maxlp
orbs2(i6)=0
enddo
!do i6=1,1
!if(orbs2(i6).eq.nn_group(i,j))goto 205
!enddo
!orbs2(1)=nn_group(i,j)
ls1=0
!print*,'k2,nn_group',k2,nn_group(i,j),i,j
if(k2.eq.nn_group(i,j))then
iii=iii+1
ls1=ls1+1
connectivity(iii,1)=k1
connectivity(iii,2)=k2
nloop(iii)=ls1+1
!print*,'nloopnloopnloop',nloop(iii)
goto 205
else
if(k1.ne.nn_group(i,j))then
!do i6=1,ii
!if(nn_group(i,j).eq.orbs1(i6))goto 205
!enddo
l=1
ii=ii+1
orbs3(ii)=nn_group(i,j)
!print*,'orbs3',orbs3(ii)
endif
endif
205 enddo

!if(l.eq.1)ls1=ls1+1
ii0=ii
!!! loop 1
do j1=1,ii0
ii=1
!print*,'loop1'
do i6=1,maxlp
orbs2(i6)=0
enddo
!do i6=1,1
!if(orbs2(i6).eq.orbs3(j1))goto 204
!enddo
orbs2(1)=orbs3(j1)
!print*,'orbs2**',(orbs2(i6),i6=1,1)
l=0
ls2=ls1
do i=1,natom
if(orbs3(j1).eq.tot_orb(i))goto 102
enddo
goto 204
102 do i1=1,nelimt(i)
!print*,'k2,nn_group',k2,nn_group(i,i1),i,i1
if(k2.eq.nn_group(i,i1))then
iii=iii+1
ls2=ls2+1
connectivity(iii,1)=k1
connectivity(iii,2)=orbs2(1)
connectivity(iii,3)=k2
nloop(iii)=ls2+1
!print*,'nloopnloopnloop',nloop(iii)
goto 204
else
if(k1.ne.nn_group(i,i1))then
!do i6=1,ii
!if(nn_group(i,i1).eq.orbs1(i6))goto 301
!enddo
l=1
ii=ii+1
orbs1(ii)=nn_group(i,i1)
!print*,'orbs1',orbs1(ii)
endif
endif

301 enddo
if(l.eq.1)ls2=ls2+l
ii1=ii
!!! loop 2
do j2=2,ii1
ii=ii1
!print*,'loop2'
do i6=2,maxlp
orbs2(i6)=0
enddo
do i6=1,1
if(orbs2(i6).eq.orbs1(j2))goto 203
enddo
orbs2(2)=orbs1(j2)
!print*,'orbs2',(orbs2(i6),i6=1,2)
l=0
ls3=ls2
do i=1,natom
if(orbs1(j2).eq.tot_orb(i))goto 103
enddo

goto 204
103 do i1=1,nelimt(i)
!print*,'k2,nn_group',k2,nn_group(i,i1),i,i1
if(k2.eq.nn_group(i,i1))then
iii=iii+1
ls3=ls3+1
connectivity(iii,1)=k1
do i6=1,2
connectivity(iii,i6+1)=orbs2(i6)
enddo
connectivity(iii,4)=k2
nloop(iii)=ls3+1
!print*,'nloopnloopnloop',nloop(iii)
goto 203
else
if(k1.ne.nn_group(i,i1))then
!do i6=1,ii
!if(nn_group(i,i1).eq.orbs1(i6))goto 302
!enddo
ii=ii+1
l=1
orbs1(ii)=nn_group(i,i1)
!print*,'orbs1',orbs1(ii)
endif
endif

302 enddo
ii2=ii
if(l.eq.1)ls3=ls3+l
!!! loop 3
do j3=ii1+1,ii2
ii=ii2
!print*,'loop3'
do i6=3,maxlp
orbs2(i6)=0
enddo
do i6=1,2
if(orbs2(i6).eq.orbs1(j3))goto 202
enddo
orbs2(3)=orbs1(j3)
!print*,'orbs2',(orbs2(i6),i6=1,3)
l=0
ls4=ls3
do i=1,natom
if(orbs1(j3).eq.tot_orb(i))goto 104
enddo

goto 203
104 do i1=1,nelimt(i)
!print*,'k2,nn_group',k2,nn_group(i,i1),i,i1
if(k2.eq.nn_group(i,i1))then
iii=iii+1
ls4=ls4+1
connectivity(iii,1)=k1
do i6=1,3
connectivity(iii,i6+1)=orbs2(i6)
enddo
connectivity(iii,5)=k2
nloop(iii)=ls4+1
!print*,'nloopnloopnloop',nloop(iii)
goto 202
else
if(k1.ne.nn_group(i,i1))then
!do i6=1,ii
!if(nn_group(i,i1).eq.orbs1(i6))goto 303
!enddo
ii=ii+1
l=1
orbs1(ii)=nn_group(i,i1)
!print*,'orbs1',orbs1(ii)
endif
endif

303 enddo
ii3=ii
if(l.eq.1)ls4=ls4+l
!!! loop 4
do j4=ii2+1,ii3
ii=ii3
!print*,'loop4'
do i6=4,maxlp
orbs2(i6)=0
enddo
do i6=1,3
if(orbs2(i6).eq.orbs1(j4))goto 201
enddo
orbs2(4)=orbs1(j4)
!print*,'orbs2**',(orbs2(i6),i6=1,4)
l=0
ls5=ls4
do i=1,natom
if(orbs1(j4).eq.tot_orb(i))goto 105
enddo
goto 202

105 do i1=1,nelimt(i)
!print*,'k2,nn_group',k2,nn_group(i,i1),i,i1
if(k2.eq.nn_group(i,i1))then
iii=iii+1
ls5=ls5+1
connectivity(iii,1)=k1
do i6=1,4
connectivity(iii,i6+1)=orbs2(i6)
enddo
connectivity(iii,6)=k2
nloop(iii)=ls5+1
!print*,'nloopnloopnloop',nloop(iii)
goto 201
else
if(k1.ne.nn_group(i,i1))then
!do i6=1,ii
!if(nn_group(i,i1).eq.orbs1(i6))goto 304
!enddo
ii=ii+1
l=1
orbs1(ii)=nn_group(i,i1)
!print*,'orbs1',orbs1(ii)
endif
endif

304 enddo
if(l.eq.1)ls5=ls5+l
ii4=ii
!!! loop 5
do j5=ii3+1,ii4
ii=ii4
!print*,'loop5',j5,ii3+1,ii4,'|',orbs1(j5),ii
do i6=5,maxlp
orbs2(i6)=0
enddo
do i6=1,4
if(orbs2(i6).eq.orbs1(j5))goto 200
enddo
orbs2(5)=orbs1(j5)
!print*,'orbs2**',(orbs2(i6),i6=1,5)
l=0
ls6=ls5
do i=1,natom
if(orbs1(j5).eq.tot_orb(i))goto 106
enddo
goto 201

106 do i1=1,nelimt(i)
!!print*,'k2,nn_group',k2,nn_group(i,i1),i,i1
if(k2.eq.nn_group(i,i1))then
iii=iii+1
ls6=ls6+1
connectivity(iii,1)=k1
do i6=1,5
connectivity(iii,i6+1)=orbs2(i6)
enddo
connectivity(iii,7)=k2
nloop(iii)=ls6+1
!print*,'nloopnloopnloop',nloop(iii)
goto 200
else
if(k1.ne.nn_group(i,i1))then
!do i6=1,ii
!if(nn_group(i,i1).eq.orbs1(i6))goto 305
!enddo
ii=ii+1
l=1
orbs1(ii)=nn_group(i,i1)
!print*,'orbs1',orbs1(ii)
endif
endif

305 enddo
if(l.eq.1)ls6=ls6+l
ii5=ii
!!! loop 6
do j6=ii4+1,ii5
ii=ii5
!print*,'loop6',ii
do i6=6,maxlp
orbs2(i6)=0
enddo
do i6=1,5
if(orbs2(i6).eq.orbs1(j6))goto 199
enddo
orbs2(6)=orbs1(j6)
!print*,'orbs2**',(orbs2(i6),i6=1,6)
l=0
ls7=ls6
do i=1,natom
if(orbs1(j6).eq.tot_orb(i))goto 107
enddo
goto 200

107 do i1=1,nelimt(i)
!!print*,'k2,nn_group',k2,nn_group(i,i1),i,i1
if(k2.eq.nn_group(i,i1))then
iii=iii+1
ls7=ls7+1
connectivity(iii,1)=k1
do i6=1,6
connectivity(iii,i6+1)=orbs2(i6)
enddo
connectivity(iii,8)=k2
nloop(iii)=ls7+1
!print*,'nloopnloopnloop',nloop(iii)
goto 199
else
if(k1.ne.nn_group(i,i1))then
!do i6=1,ii
!if(nn_group(i,i1).eq.orbs1(i6))goto 306
!enddo
ii=ii+1
l=1
orbs1(ii)=nn_group(i,i1)
!print*,'orbs1',orbs1(ii)
endif
endif

306 enddo
if(l.eq.1)ls7=ls7+l
ii6=ii
!!! loop 7
do j7=ii5+1,ii6
ii=ii6
!print*,'loop7',ii
do i6=7,maxlp
orbs2(i6)=0
enddo
do i6=1,6
if(orbs2(i6).eq.orbs1(j7))goto 198
enddo
orbs2(7)=orbs1(j7)
!print*,'orbs2**',(orbs2(i6),i6=1,7)
l=0
ls8=ls7
do i=1,natom
if(orbs1(j7).eq.tot_orb(i))goto 108
enddo
goto 199

108 do i1=1,nelimt(i)
!!print*,'k2,nn_group',k2,nn_group(i,i1),i,i1
if(k2.eq.nn_group(i,i1))then
iii=iii+1
ls8=ls8+1
connectivity(iii,1)=k1
do i6=1,7
connectivity(iii,i6+1)=orbs2(i6)
enddo
connectivity(iii,9)=k2
nloop(iii)=ls8+1
!print*,'nloopnloopnloop',nloop(iii)
goto 198
else
if(k1.ne.nn_group(i,i1))then
!do i6=1,ii
!if(nn_group(i,i1).eq.orbs1(i6))goto 307
!enddo
ii=ii+1
l=1
orbs1(ii)=nn_group(i,i1)
!print*,'orbs1',orbs1(ii)
endif
endif

307 enddo
if(l.eq.1)ls8=ls8+l
ii7=ii
!!! loop 8
do j8=ii6+1,ii7
ii=ii7
!print*,'loop8',ii
do i6=8,maxlp
orbs2(i6)=0
enddo
do i6=1,7
if(orbs2(i6).eq.orbs1(j8))goto 197
enddo
orbs2(8)=orbs1(j8)
!print*,'orbs2**',(orbs2(i6),i6=1,8)
l=0
ls9=ls8
do i=1,natom
if(orbs1(j8).eq.tot_orb(i))goto 109
enddo
goto 198

109 do i1=1,nelimt(i)
!print*,'k2,nn_group',k2,nn_group(i,i1),i,i1
if(k2.eq.nn_group(i,i1))then
iii=iii+1
ls9=ls9+1
connectivity(iii,1)=k1
do i6=1,8
connectivity(iii,i6+1)=orbs2(i6)
enddo
connectivity(iii,10)=k2
nloop(iii)=ls9+1
!print*,'nloopnloopnloop',nloop(iii)
goto 197
else
if(k1.ne.nn_group(i,i1))then
!do i6=1,ii
!if(nn_group(i,i1).eq.orbs1(i6))goto 308
!enddo
ii=ii+1
l=1
orbs1(ii)=nn_group(i,i1)
!print*,'orbs1',orbs1(ii)
endif
endif

308 enddo
if(l.eq.1)ls9=ls9+l
ii8=ii
!!! loop 9
do j9=ii7+1,ii8
ii=ii8
!print*,'loop9',ii
do i6=9,maxlp
orbs2(i6)=0
enddo
do i6=1,8
if(orbs2(i6).eq.orbs1(j9))goto 196
enddo
orbs2(9)=orbs1(j9)
!print*,'orbs2**',(orbs2(i6),i6=1,9)
l=0
ls10=ls9
do i=1,natom
if(orbs1(j9).eq.tot_orb(i))goto 110
enddo
goto 197

110 do i1=1,nelimt(i)
!print*,'k2,nn_group',k2,nn_group(i,i1),i,i1
if(k2.eq.nn_group(i,i1))then
iii=iii+1
ls10=ls10+1
connectivity(iii,1)=k1
do i6=1,9
connectivity(iii,i6+1)=orbs2(i6)
enddo
connectivity(iii,11)=k2
nloop(iii)=ls10+1
!print*,'nloopnloopnloop',nloop(iii)
goto 196
else
if(k1.ne.nn_group(i,i1))then
!do i6=1,ii
!if(nn_group(i,i1).eq.orbs1(i6))goto 309
!enddo
ii=ii+1
l=1
orbs1(ii)=nn_group(i,i1)
!print*,'orbs1',orbs1(ii)
endif
endif

309 enddo
if(l.eq.1)ls10=ls10+l
ii9=ii
!!! loop 10
do j10=ii8+1,ii9
ii=ii9
!print*,'loop10'
do i6=10,maxlp
orbs2(i6)=0
enddo
do i6=1,9
if(orbs2(i6).eq.orbs1(j10))goto 195
enddo
orbs2(10)=orbs1(j10)
!print*,'orbs2**',(orbs2(i6),i6=1,10)
l=0
ls11=ls10
do i=1,natom
if(orbs1(j10).eq.tot_orb(i))goto 111
enddo
goto 196

111 do i1=1,nelimt(i)
!print*,'k2,nn_group',k2,nn_group(i,i1),i,i1
if(k2.eq.nn_group(i,i1))then
iii=iii+1
ls11=ls11+1
connectivity(iii,1)=k1
do i6=1,10
connectivity(iii,i6+1)=orbs2(i6)
enddo
connectivity(iii,12)=k2
nloop(iii)=ls11+1
!print*,'nloopnloopnloop',nloop(iii)
goto 195
else
if(k1.ne.nn_group(i,i1))then
!do i6=1,ii
!if(nn_group(i,i1).eq.orbs1(i6))goto 310
!enddo
ii=ii+1
l=1
orbs1(ii)=nn_group(i,i1)
!print*,'orbs1',orbs1(ii)
endif
endif

310 enddo
if(l.eq.1)ls11=ls11+l
ii10=ii
!!! loop 11
do j11=ii9+1,ii10
ii=ii10
!print*,'loop11'
do i6=11,maxlp
orbs2(i6)=0
enddo
do i6=1,10
if(orbs2(i6).eq.orbs1(j11))goto 194
enddo
orbs2(11)=orbs1(j11)
!print*,'orbs2**',(orbs2(i6),i6=1,11)
l=0
ls12=ls11
do i=1,natom
if(orbs1(j11).eq.tot_orb(i))goto 112
enddo
goto 195

112 do i1=1,nelimt(i)
!print*,'k2,nn_group',k2,nn_group(i,i1),i,i1
if(k2.eq.nn_group(i,i1))then
iii=iii+1
ls12=ls12+1
connectivity(iii,1)=k1
do i6=1,11
connectivity(iii,i6+1)=orbs2(i6)
enddo
connectivity(iii,13)=k2
nloop(iii)=ls12+1
!print*,'nloopnloopnloop',nloop(iii)
goto 194
else
if(k1.ne.nn_group(i,i1))then
!do i6=1,ii
!if(nn_group(i,i1).eq.orbs1(i6))goto 311
!enddo
ii=ii+1
l=1
orbs1(ii)=nn_group(i,i1)
!print*,'orbs1',orbs1(ii),ls1
endif
endif

311 enddo
if(l.eq.1)ls12=ls12+l
ii11=ii
!!! loop 12
!print*,'loop12'
do j12=ii10+1,ii11
ii=ii11
do i6=12,maxlp
orbs2(i6)=0
enddo
do i6=1,11
if(orbs2(i6).eq.orbs1(j12))goto 193
enddo
orbs2(12)=orbs1(j12)
!print*,'orbs2**',(orbs2(i6),i6=1,12)
l=0
ls13=ls12
do i=1,natom
if(orbs1(j12).eq.tot_orb(i))goto 113
enddo
goto 194

113 do i1=1,nelimt(i)
!!print*,'k2,nn_group',k2,nn_group(i,i1),i,i1
if(k2.eq.nn_group(i,i1))then
iii=iii+1
ls13=ls13+1
connectivity(iii,1)=k1
do i6=1,12
connectivity(iii,i6+1)=orbs2(i6)
enddo
connectivity(iii,14)=k2
nloop(iii)=ls13+1
!print*,'nloopnloopnloop',nloop(iii)
goto 193
else
if(k1.ne.nn_group(i,i1))then
!do i6=1,ii
!if(nn_group(i,i1).eq.orbs1(i6))goto 312
!enddo
ii=ii+1
l=1
orbs1(ii)=nn_group(i,i1)
!print*,'orbs1',orbs1(ii)
endif
endif

312 enddo
if(l.eq.1)ls13=ls13+l
ii12=ii
!!! loop 13
!print*,'loop13'
do j13=ii11+1,ii12
ii=ii12
do i6=13,maxlp
orbs2(i6)=0
enddo
do i6=1,12
if(orbs2(i6).eq.orbs1(j13))goto 192
enddo
orbs2(13)=orbs1(j13)
!print*,'orbs2**',(orbs2(i6),i6=1,13)
l=0
ls14=ls13
do i=1,natom
if(orbs1(j13).eq.tot_orb(i))goto 114
enddo
goto 193

114 do i1=1,nelimt(i)
!print*,'k2,nn_group',k2,nn_group(i,i1),i,i1
if(k2.eq.nn_group(i,i1))then
iii=iii+1
ls14=ls14+1
connectivity(iii,1)=k1
do i6=1,13
connectivity(iii,i6+1)=orbs2(i6)
enddo
connectivity(iii,15)=k2
nloop(iii)=ls14+1
!print*,'nloopnloopnloop',nloop(iii)
goto 192
else
if(k1.ne.nn_group(i,i1))then
!do i6=1,ii
!if(nn_group(i,i1).eq.orbs1(i6))goto 313
!enddo
ii=ii+1
orbs1(ii)=nn_group(i,i1)
!print*,'orbs1',orbs1(ii)
endif
endif

313 enddo
192 enddo

193 enddo

194 enddo

195 enddo

196 enddo

197 enddo

198 enddo

199 enddo

200 enddo

201 enddo

202 enddo

203 enddo

204 enddo

!205 enddo

n=1000
do i1=1,iii
if(n.gt.nloop(i1))n=nloop(i1)
enddo

loop_score=n
nn=0
l=0
do i1=1,iii
if(loop_score.eq.nloop(i1))then
l=l+1
nnn(l)=i1
nn=nn+1
endif
enddo

print*,'nnn**',(nnn(i),i=1,l)

do i2=1,20
do i1=1,20
connectivity_row(i1,i2)=0
enddo
enddo
!print*,'',nnn,nloop(nnn)
!if (nn.eq.1)then
do i2=1,l
do i1=1,nloop(nnn(i2))+1
connectivity_row(i2,i1)=connectivity(nnn(i2),i1)
!print*,'connectivity_row',connectivity_row(i1)
enddo
enddo

nrow=l
!endif

print*,'nloop',(nloop(i1),i1=1,iii)
do i1=1,iii
print*,'connectivity',(connectivity(i1,i2),i2=1,nloop(i1)+1)
enddo
do i2=1,l
print*,'connectivity_row',(connectivity_row(i2,i1),i1=1,nloop(nnn(i2))+1)
enddo

print*,'shortestgest loop',n
!if(k1.eq.24.and.k2.eq.29)stop

100 return
end subroutine symm_loops
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nnat_bond_cal(nl,str1,ncqs,bondq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! This subroutine calculates the number of nearest neighbour bonds !!!!!

use commondat
implicit none

integer::m19,m18,i1,i2,i3,i4,i5,i6,i7,ii,iii,iiii,nl,ncqs
integer::nn(10),str1(15000,20),bondq(15000)

print*,'enter nnat_bond_cal'
do i1=1,15000
bondq(i1)=0
enddo

if(nnnatom.eq.0) goto 900
231 format (30I3)
i4=0
do i1=1,ncqs
write(*,231)(str1(i1,i2),i2=1,nae)
iiii=1+(nae-nl*2-nlast)/2
do i3=1+nl*2,(nae-nlast),2
iii=0
nn(1)=0
nn(2)=0
do i5=i3,i3+1
do i4=1,atom
do i7=1,atn(active_atoms(i4))
if(str1(i1,i5).eq.atoset(active_atoms(i4),i7))then
iii=iii+1
if(mod(i5,2).eq.0) nn(2)=i4
if(mod(i5,2).eq.1) nn(1)=i4
goto 517
endif
enddo

enddo
517 enddo
if(nn(2).eq.nn(1))goto 400
if(iii.ne.2) goto 400
do i6=1,nnnatom
ii=0
do i7=1,2
do i4=1,2
if(nn(i4).eq.nnat_bond(i6,i7))then
ii=ii+1
endif
enddo
enddo
if(ii.eq.2)then
iiii=iiii-1
goto 400
endif
enddo

400 enddo
bondq(i1)=iiii
!print*,'nnbd:iiii',iiii
enddo

print*,'exit nnat_bond_cal'
900 return
end subroutine nnat_bond_cal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nnat_bond_cal_2(nl,str1,ncqs,bondq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::m19,m18,i,i1,i2,i3,i4,i5,i6,i7,iii,iiii,nl,ncqs,k
integer::nn(10),str1(15000,20),bondq(15000),sl(5000)
real::least,bondq1_dist(15000)
double precision::ii,bondq_dist(15000)

print*,'enter nnat_bond_cal_2'

open(unit=13,file='nnbd.temp',status='unknown')
do i1=1,15000
bondq(i1)=0
bondq_dist(i1)=0.0
enddo

231 format (30I3)
i4=0
do i1=1,ncqs
ii=0.0
!write(*,231)i1,(str1(i1,i2),i2=1,nae)
do i3=1+nl*2,(nae-nlast),2
iii=0
nn(1)=0
nn(2)=0
do i5=i3,i3+1
do i4=1,atom
do i7=1,atn(active_atoms(i4))
!print*,'atn(i4)',atn(active_atoms(i4)),active_atoms(i4)
!print*,'atoset',str1(i1,i5),atoset(active_atoms(i4),i7)
if(str1(i1,i5).eq.atoset(active_atoms(i4),i7))then
iii=iii+1
if(mod(i5,2).eq.0) nn(2)=active_atoms(i4)
if(mod(i5,2).eq.1) nn(1)=active_atoms(i4)
!print*,'------------',nn(1),nn(2)
goto 517
endif
enddo

enddo
517 enddo
print*,'nn(2),nn(1)',nn(2),nn(1),iii
if(iii.ne.2) goto 400
if(nn(2).eq.nn(1))then

!!! changes done for EDEN
iab_length=1.00
!!!!!!!!!!!!!!!!!!!!!!!!!

if (iab_length.eq.0.0)ii=ii+100.0
if (iab_length.ne.0.0)ii=ii+iab_length
!print*,'iab_length',iab_length
goto 400
endif
ii=ii+dist_act_rel_mat(nn(1),nn(2))
!print*,'dist_act_rel_mat(nn(1),nn(2))',dist_act_rel_mat(nn(1),nn(2)),ii
400 enddo
bondq_dist(i1)=ii
!print*,'bondq_dist(i1)',bondq_dist(i1),ii
write(13,*)bondq_dist(i1)
enddo
rewind(13)
do i1=1,ncqs
read(13,105)bondq1_dist(i1)
!write(*,*),i1,bondq1_dist(i1)
enddo
105 format (F8.4)
iii=0
iiii=0
211 least=10000.0
do i=1,ncqs
do i1=1,iii
if(i.eq.sl(i1))goto 210
enddo
if(bondq1_dist(i).le.least)least=bondq1_dist(i)
210 enddo
!print*,'least',least
if(least.ne.100.0.or.least.ne.0.0)iiii=iiii+1
do i=1,ncqs
if(bondq1_dist(i).eq.least)then
iii=iii+1
bondq(i)=iiii
sl(iii)=i
endif
enddo
if(iii.ne.ncqs)goto 211

!do i=1,ncqs
!write(*,231)i,(str1(i,i2),i2=1,nae),bondq(i)
!!print*,'bondq(i)',i,bondq(i)
!enddo
CALL SYSTEM ("rm nnbd.temp")
print*,'exit nnat_bond_cal_2'
900 return
end subroutine nnat_bond_cal_2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine quality_factor(nl,str1,ncqs,q_fac,str_quality_1,str_quality_2,bondq)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/ats/atsymset,nsym,syn,at_sym

integer::i,j,nl,qt,qt_cnt(10),ncqs,i9,m16,m17,m18,m19,qtg,str1(15000,20),equality,noqt&
,jj,nqul1,nqul2,i4,q_fac(15000),k6,k7,k8,str_quality_1(15000),str_quality_2(15000)&
,bondq(15000),qual_mat(6,15000),q,pref_radical(15000),qt_max,mbondq(15000)
integer::atsymset(20,20),nsym,syn(50),at_sym(50),fchrgq(15000)


print*,'itb,syb,nnb,radical,mnbond',itb,syb,nnb,radical,mnbond
if(qflg.eq.1)then 
do j=1,ncqs
q_fac(j)=1
enddo
goto 555
endif
if(itb.ne.1.or.flg1.eq.1) call intra_bond_factor(nl,str1,ncqs,str_quality_1)
if(syb.ne.1.and.nsym.ne.1) call symm_break_factor(nl,str1,ncqs,str_quality_2)
if(input_flg.eq.1)then
  if(nnb.ne.1.or.flg1.eq.1) call nnat_bond_cal(nl,str1,ncqs,bondq)
endif
if(input_flg.eq.0)then
  if(nnb.ne.1.or.flg1.eq.1) call nnat_bond_cal_2(nl,str1,ncqs,bondq)
endif
if(radical.ne.1.and.prad.ne.0.and.nlast.ne.0.or.flg1.eq.1) call prio_rad_str(nl,str1,ncqs,pref_radical)
if(mnbond.ne.1.and.imbd.ne.0.or.flg1.eq.1) call main_bond_cal(nl,str1,ncqs,mbondq)

do j=1,ncqs
  do i=1,6
  qual_mat(i,j)=0
  enddo
enddo

do j=1,ncqs
i=itb
qual_mat(i,j)=qual_mat(i,j)+str_quality_1(j)
i=syb
qual_mat(i,j)=qual_mat(i,j)+str_quality_2(j)
i=nnb
qual_mat(i,j)=qual_mat(i,j)+bondq(j)
i=radical
qual_mat(i,j)=qual_mat(i,j)+pref_radical(j)
i=mnbond
qual_mat(i,j)=qual_mat(i,j)+mbondq(j)
enddo

213 format(10I5)

noqt=1
do i=2,6
  if(qual_mat(i,1).gt.0)then
  noqt=noqt+1
  endif
enddo

if(noqt.gt.2)then
  do i=2,noqt-1
  qt_max=1
    do j=1,ncqs
      if(qt_max.lt.qual_mat(i+1,j))qt_max=qual_mat(i+1,j)
    enddo
    do j=1,ncqs
      qual_mat(i+1,j)=qt_max*(qual_mat(i,j)-1)+qual_mat(i+1,j)
    enddo
  enddo
endif

do j=1,ncqs
q_fac(j)=qual_mat(noqt,j)
enddo
!do j=1,ncqs
!write(*,200),j,q_fac(j),(str1(j,i),i=1,nae)
!enddo
!stop
200 format(20I3)
555 return
end subroutine quality_factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vector_rep(nl,str1,totstr,vec)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! vector representation of the structures for independency calculations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none 

common/ats/atsymset,nsym,syn,at_sym
integer::j,j1,j2,i,i1,i2,i3,i4,i5,i6,i7,i8,k,k6,k7,k8,k11,k12,k13,k14,k15,l,nend,c1,c2,c3,c4,c5,c6,c,d,e,f,totstr,&
cof,ntdet,det,det1,nl,permtn1,permtn2,permtn,tndet,kk,kkk,x,y,z,detno,n1,n2,n3,n4,n5,n6,n7,n8,i9,i10,i11,&
c7,c8,c9,c10
integer::rum1(100),rum2(100),m(500),bonds(100),coef(500),vec(15000,1000),str1(15000,20)
integer::score(40),score_sym(40),col(20),detcount,col1(2,20)
real*8::factorial,norm_const,detn
Double Precision::D1(1000),aldet_new,aldet(1000)
integer::atsymset(20,20),nsym,syn(50),at_sym(50)


   integer,dimension(:,:),allocatable::str2,beta,beta1,alfa,vector,all_det,sign_det,all_det_new

print*,'enter vector_rep'
x=100000
y=10000
z=1000

allocate(sign_det(x,z))
allocate(str2(x,y))
allocate(beta(x,y))
allocate(beta1(x,y))
allocate(alfa(x,y))
allocate(vector(x,y))
allocate(all_det(x,y))
allocate(all_det_new(x,y))

str2(x,y)=0
beta(x,y)=0
beta1(x,y)=0
alfa(x,y)=0
vector(x,y)=0
all_det(x,y)=0
str_det_sec(x,z)=0

do i=1,10000
strdet(i)=0
enddo

print*,'sub: vector 1'
!!!!!!! "rum1()" is the orbital numbers except lone paires !!!!!
!!!!!!! "rum2()" is the orbital numbers except lone paires arranged in increasing order !!!!!

!do i1=1,totstr
!write(10,*),(str1(i1,i2),i2=1,nae)
!enddo

do k6=1,100
rum1(k6)=0
enddo
do k6=1,nae-nl*2
!print*,'nnlp',k6
do k7=nl*2+1,nae
do k8=1,k6
if(rum1(k8).eq.str1(1,k7))goto 469
enddo
if(str1(1,k7).gt.rum1(k6))then
rum1(k6)=str1(1,k7)
!print*,rum1(k6)
endif
469 enddo
enddo
do k6=1,100
rum2(k6)=0
enddo
k7=0
do k6=nae-nl*2,1,-1
k7=k7+1
rum2(k7)=rum1(k6)

enddo
nend=k7+mod(k7,2)


!!!!!!! "all_det()" are the all possible determinants of the set. Ex. for 5
!electrons in 5 orbitals case : there must be 3 alfa and two beta in each
!determinant so the possible number of determinants are= !5/(!3*!2)=10 so all_det()
!is the 10 defferent sets of determinant. 

k11=0
k12=0
do i1=1,k7
!print*,'sourav'
k11=1
m(1)=i1
if (k11.eq.nend/2) then
k12=k12+1
do k15=1,nend/2
all_det(k12,k15)=rum2(m(k15))
enddo
goto 100
endif

do i2=i1,k7
k11=2
if(i2.eq.m(1))goto 101
m(2)=i2
if (k11.eq.nend/2) then
k12=k12+1
do k15=1,nend/2
all_det(k12,k15)=rum2(m(k15))
enddo
goto 101
endif

do i3=i2,k7
k11=3
do k13=1,2
if(i3.eq.m(k13))goto 102
enddo
m(3)=i3
if (k11.eq.nend/2) then
k12=k12+1
do k15=1,nend/2
all_det(k12,k15)=rum2(m(k15))
enddo
goto 102
endif

do i4=i3,k7
k11=4
do k14=1,3
if(i4.eq.m(k14))goto 103
enddo
m(4)=i4
if (k11.eq.nend/2) then
k12=k12+1
do k15=1,nend/2
all_det(k12,k15)=rum2(m(k15))
enddo
goto 103
endif

!print*,'sourav4'
do i5=i4,k7
k11=5
do k14=1,4
if(i5.eq.m(k14))goto 104
enddo
m(5)=i5
if (k11.eq.nend/2) then
k12=k12+1
do k15=1,nend/2
all_det(k12,k15)=rum2(m(k15))
enddo
goto 104
endif

do i6=i5,k7
k11=6
do k14=1,5
if(i6.eq.m(k14))goto 105
enddo
m(6)=i6
if (k11.eq.nend/2) then
k12=k12+1
do k15=1,nend/2
all_det(k12,k15)=rum2(m(k15))
enddo
!print*,(all_det(k12,l),l=1,k7/2)
goto 105
endif

do i7=i6,k7
k11=7
do k14=1,6
if(i7.eq.m(k14))goto 106
enddo
m(7)=i7
if (k11.eq.nend/2) then
k12=k12+1
do k15=1,nend/2
all_det(k12,k15)=rum2(m(k15))
enddo
!print*,(all_det(k12,l),l=1,k7/2)
goto 106
endif

do i8=i7,k7
k11=8
do k14=1,7
if(i8.eq.m(k14))goto 107
enddo
m(8)=i8
if (k11.eq.nend/2) then
k12=k12+1
do k15=1,nend/2
all_det(k12,k15)=rum2(m(k15))
enddo
!print*,(all_det(k12,l),l=1,k7/2)
goto 107
endif

do i9=i8,k7
k11=9
do k14=1,8
if(i9.eq.m(k14))goto 108
enddo
m(9)=i9
if (k11.eq.nend/2) then
k12=k12+1
do k15=1,nend/2
all_det(k12,k15)=rum2(m(k15))
enddo
!print*,(all_det(k12,l),l=1,k7/2)
goto 108
endif

do i10=i9,k7
k11=10
do k14=1,9
if(i10.eq.m(k14))goto 109
enddo
m(10)=i10
if (k11.eq.nend/2) then
k12=k12+1
do k15=1,nend/2
all_det(k12,k15)=rum2(m(k15))
enddo
!print*,(all_det(k12,l),l=1,k7/2)
goto 109
endif

109 enddo
108 enddo
107 enddo
106 enddo
105 enddo
104 enddo
103 enddo
102 enddo
101 enddo
100 enddo

ntdet=k12
!print*,'sub: vector 2'

d=nae-nl*2
e=nend/2
f=d-e
c=factorial(d)/(factorial(e)*factorial(f))
do i1=1,c
i4=nend/2
do i3=1,nae-nl*2
do i2=1,nend/2
if(all_det(i1,i2).eq.rum2(i3))goto 557
enddo
i4=i4+1
!print*,rum2(i3),i4
all_det(i1,i4)=rum2(i3)
557 enddo
nbeta=(nae-(nl*2)-(mult-1))/2
nalpha=nbeta+(mult-1)

if (nl.ne.0)then
i2=0
do i3=1,nalpha
i2=i2+1
all_det_new(i1,i2)=all_det(i1,i3)

enddo
do i3=1,nl,2
i2=i2+1
all_det_new(i1,i2)=str1(1,i3)
enddo
do i3=nalpha+1,nalpha+nbeta
i2=i2+1
all_det_new(i1,i2)=all_det(i1,i3)
enddo
do i3=1,nl,2
i2=i2+1
all_det_new(i1,i2)=str1(1,i3)
!print*,'all_det_new(i1,i4)=str1(1,i3)',i2,all_det_new(i1,i2),str1(1,i3)
enddo
endif

if(nl.eq.0)then
do i3=1,nae
all_det_new(i1,i3)=all_det(i1,i3)
enddo
endif

!print*,i1,'all_det',(all_det(i1,i2),i2=1,nae-nl*2)
!write(9,300),'all_det_new',i1,(all_det_new(i1,i2),i2=1,nalpha+nbeta+nl*2)

enddo
tot_ndet=i1
do i1=1,15000
do i2=1,1000
vector(i1,i2)=0
enddo
enddo

i4=0
tndet=0
do i7=1,totstr
!write(9,*),(str1(i7,i2),i2=1,nae)
j=0
i2=0
detcount=0
do i4=1,100
bonds(100)=0
enddo
do i1=nl*2+1,nae-nlast
i2=i2+1
bonds(i2)=str1(i7,i1)
enddo

!!!!!!  'struc_albe () ' is the determinants of each structures presenting only alfas !!!!!!
!!!!!!  for two bonds and one unpaired electrons Ex. (1-2)(3-4)5, two bonds
!!!!!!  between orb 1,2 and orbs 3,4 and unpaired electron in in 5. then the
!!!!!!  determinants are (a1b2-a2b1)(a3b4-a4b3)a5 taking unpaired as alfa. so the
!!!!!!  determinants are a1b2a3b4a5-a1b2a4b3a5-a2b1a3b4a5+a2b1a4b3a5 to delet the a,b
!!!!!!  notation or to make it in only number sets we put all alfas first and beta
!!!!!!  second, therefore the number set representing the above determinants are =
!!!!!!  135/24 - 145/23 - 235/14 + 245/13, 'cof()' is the sign of the determinants.  !!!!!!!

do i1=1,500
m(i1)=0
enddo

print*,'sub: vector 3'
k11=0
k12=0
cof=1
do i1=1,2
cof=cof*(-1)**(i1-1)
c1=cof
k11=1
m(1)=i1
if (k11.eq.(nae-nl*2-nlast)/2) then
k12=k12+1
do k15=1,(nae-nl*2-nlast)/2
str2(k12,k15)=bonds(m(k15))
coef(k12)=cof
enddo
print*,'sourav1',coef(k12),(str2(k12,l),l=1,(nae-nl*2-nlast)/2)
goto 200
endif

do i2=3,4
cof=cof*(-1)**(i2-1)
c2=cof
k11=2
if(i2.eq.m(1))goto 201
m(2)=i2
if (k11.eq.(nae-nl*2-nlast)/2) then
k12=k12+1
do k15=1,(nae-nl*2-nlast)/2
str2(k12,k15)=bonds(m(k15))
coef(k12)=cof
enddo
print*,'sourav2',coef(k12),(str2(k12,l),l=1,(nae-nl*2-nlast)/2)
goto 201
endif

do i3=5,6
cof=cof*(-1)**(i3-1)
c3=cof
k11=3
!print*,'sourav1'
do k13=1,2
if(i3.eq.m(k13))goto 202
enddo
!print*,'sourav2'
m(3)=i3
if (k11.eq.(nae-nl*2-nlast)/2) then
k12=k12+1
do k15=1,(nae-nl*2-nlast)/2
str2(k12,k15)=bonds(m(k15))
coef(k12)=cof
enddo
!print*,'sourav',coef(k12),(str2(k12,l),l=1,(nae-nl*2-nlast)/2)
goto 202
endif

do i4=7,8
cof=cof*(-1)**(i4-1)
c4=cof
k11=4
do k14=1,3
if(i4.eq.m(k14))goto 203
enddo
m(4)=i4
!print*,(m(l),l=1,nao)
if (k11.eq.(nae-nl*2-nlast)/2) then
k12=k12+1
do k15=1,(nae-nl*2-nlast)/2
str2(k12,k15)=bonds(m(k15))
coef(k12)=cof
enddo
!print*,'sourav',coef(k12),(str2(k12,l),l=1,(nae-nl*2-nlast)/2)
goto 203
endif

!print*,'sourav4'
do i5=9,10
cof=cof*(-1)**(i5-1)
c5=cof
k11=5
do k14=1,4
if(i5.eq.m(k14))goto 204
enddo
m(5)=i5
if (k11.eq.(nae-nl*2-nlast)/2) then
k12=k12+1
do k15=1,(nae-nl*2-nlast)/2
str2(k12,k15)=bonds(m(k15))
coef(k12)=cof
enddo
!print*,(strc(i,l),l=1,nae-nl*2-nlast)
goto 204
endif

do i6=11,12
cof=cof*(-1)**(i6-1)
c6=cof
k11=6
do k14=1,5
if(i6.eq.m(k14))goto 205
enddo
m(6)=i6
if (k11.eq.(nae-nl*2-nlast)/2) then
k12=k12+1
do k15=1,(nae-nl*2-nlast)/2
str2(k12,k15)=bonds(m(k15))
coef(k12)=cof
enddo
!print*,(str2(i,l),l=1,nae-nlp*2-nlast)
goto 205
endif

do i8=13,14
cof=cof*(-1)**(i8-1)
c7=cof
k11=7
do k14=1,6
if(i8.eq.m(k14))goto 206
enddo
m(7)=i8
if (k11.eq.(nae-nl*2-nlast)/2) then
k12=k12+1
do k15=1,(nae-nl*2-nlast)/2
str2(k12,k15)=bonds(m(k15))
coef(k12)=cof
enddo
!print*,(str2(i,l),l=1,nae-nlp*2-nlast)
goto 206
endif

do i9=15,16
cof=cof*(-1)**(i9-1)
c8=cof
k11=8
do k14=1,7
if(i9.eq.m(k14))goto 207
enddo
m(8)=i9
if (k11.eq.(nae-nl*2-nlast)/2) then
k12=k12+1
do k15=1,(nae-nl*2-nlast)/2
str2(k12,k15)=bonds(m(k15))
coef(k12)=cof
enddo
!print*,(str2(i,l),l=1,nae-nlp*2-nlast)
goto 207
endif

do i10=17,18
cof=cof*(-1)**(i10-1)
c9=cof
k11=9
do k14=1,8
if(i10.eq.m(k14))goto 208
enddo
m(9)=i10
if (k11.eq.(nae-nl*2-nlast)/2) then
k12=k12+1
do k15=1,(nae-nl*2-nlast)/2
str2(k12,k15)=bonds(m(k15))
coef(k12)=cof
enddo
goto 208
endif

do i11=19,20
cof=cof*(-1)**(i11-1)
c10=cof
k11=10
do k14=1,9
if(i11.eq.m(k14))goto 209
enddo
m(10)=i11
if (k11.eq.(nae-nl*2-nlast)/2) then
k12=k12+1
do k15=1,(nae-nl*2-nlast)/2
str2(k12,k15)=bonds(m(k15))
coef(k12)=cof
enddo
goto 209
endif

209 cof=c10
enddo
208 cof=c9
enddo
207 cof=c8
enddo
206 cof=c7
enddo
205 cof=c6
enddo
204 cof=c5
enddo
203 cof=c4
enddo
202 cof=c3
enddo
201 cof=c2
enddo
200 cof=c1
enddo

ndet=k12

do k14=1,k12
tndet=tndet+1
i8=(nae-nl*2-nlast)/2
do i6=nae-nlast+1,nae
i8=i8+1
str2(k14,i8)=str1(i7,i6)
enddo
i5=0
do i1=1,i8
do i2=nl*2+1,nae-nlast
if(str2(k14,i1).eq.str1(i7,i2))then
if(mod(i2,2).eq.1)i3=i2+1
if(mod(i2,2).eq.0)i3=i2-1
endif
enddo
beta1(k14,i1)=str1(i7,i3)
enddo

k7=0
if(nl.ne.0)then
do k6=1,nl*2,2
k7=k7+1
alfa(k14,k7)=str1(i7,k6)
beta(k14,k7)=str1(i7,k6)
enddo
endif

do k6=1,(nae-nl*2-nlast)/2
k7=k7+1
alfa(k14,k7)=str2(k14,k6)
beta(k14,k7)=beta1(k14,k6)
enddo

!write(9,*),'det:vector',(alfa(k14,k6),k6=1,k7),(beta(k14,k6),k6=1,k7)
!print*,'sub: vector 4'

if(nlast.ne.0)then
do k6=nae-nlast+1,nae
k7=k7+1
alfa(k14,k7)=str1(i7,k6)
enddo
endif

do k6=1,100
rum1(k6)=0
enddo

k7=0
i5=0
kkk=0
permtn1=0
770 k=100
kk=0
do k6=1,(nae-nl*2-nlast)/2+nl+nlast
do k8=1,k7
if(rum1(k8).eq.alfa(k14,k6))goto 769
enddo
kk=kk+1
if(k.gt.alfa(k14,k6))then
kkk=kk
k=alfa(k14,k6)
endif
769 enddo
k7=k7+1
i5=i5+1
rum1(k7)=k
permtn1=permtn1+(kkk-1)
detmnt(tndet,i5)=k
strdet(tndet)=i7

!print*,'strdet',strdet(tndet),tndet
!print*,'detmnt1',detmnt(tndet,i5)
!print*,'kk,k7,permtn1',kk,k7,permtn1
!print*,'rum1',rum1(k7)

if(k7.lt.((nae-nl*2-nlast)/2)+nl+nlast)goto 770


do k6=1,100
rum1(k6)=0
enddo

k7=0
permtn2=0
kkk=0
780 k=1000
kk=0
do k6=1,(nae-nl*2-nlast)/2+nl
!print*,'k6',k6
do k8=1,k7
if(rum1(k8).eq.beta(k14,k6))goto 789
enddo
kk=kk+1
!print*,'nnlp',beta(k14,k6),k
if(k.gt.beta(k14,k6))then
kkk=kk
k=beta(k14,k6)
!print*,'kk',kk
endif
789 enddo
!print*,'kkk',k
k7=k7+1
i5=i5+1
rum1(k7)=k
permtn2=permtn2+(kkk-1)
detmnt(tndet,i5)=k
strdet(tndet)=i7

!print*,'strdet',strdet(tndet),tndet
!print*,'detmnt2',detmnt(tndet,i5)

!print*,'kk,permtn2',kk,permtn2
!print*,'rum1',rum1(k7)

if(k7.lt.((nae-nl*2-nlast)/2)+nl+nlast)goto 780

permtn=permtn1+permtn2

det_sign(tndet)=(-1)**permtn
j=j+1
sign_det(i7,j)=det_sign(tndet)
detcount=detcount+1
str_det_sec(i7,detcount)=tndet

enddo
300 format (a,20I3)

!print*,'sub: vector 5'



do i1=1,k12
i4=nend/2
do i3=1,nae-nl*2
do i2=1,nend/2
if(str2(i1,i2).eq.rum2(i3))goto 558
enddo
i4=i4+1
!print*,rum2(i3),i4
str2(i1,i4)=rum2(i3)
558 enddo
enddo

i3=0
do i1=1,ndet
do i2=1,ntdet
det=0
det1=0
do i3=1,nend/2
do i4=1,nend/2
if(str2(i1,i3).eq.all_det(i2,i4))then
det=det+1
endif
enddo
enddo
if(det.eq.nend/2)then
do i3=1+nend/2,nae-nl*2
do i4=1+nend/2,nae-nl*2
if(str2(i1,i3).eq.all_det(i2,i4))then
det1=det1+1
endif
enddo
enddo
if(det1.eq.(nae-nl*2)-(nend/2))then
i3=i3+1
vector(i7,i2)=1*coef(i1)
endif
endif
enddo
enddo

i3=0
do i2=1,ntdet
if(vector(i7,i2).ne.0)then
i3=i3+1
vec(i7,i3)=i2
endif
enddo

detno=i3

detn= real(detno,4)
!write(*,321),i7,'vec',(vec(i7,i4),i4=1,i3)
norm_const=1/sqrt(detn)

!write(9,321),i7,'vector>',(vector(i7,i4),i4=1,ntdet)

enddo

!write(9,*),'ntdet,tndet,ndet',ntdet,tndet,ndet
!do i=1,totstr
!write(9,321),i,'vector>',(vector(i,i4),i4=1,ntdet)
!enddo

do l=1,totstr*ndet
!write(9,123),'detmnt(l,i)',(detmnt(l,i),i=1,nae),det_sign(l)
enddo
123 format(a,2x,50I5)
321 format(I5,a,500I3)
!!!end of independency test

339 print*,'exit vector_rep'
deallocate(sign_det)
deallocate(str2)
deallocate(beta)
deallocate(beta1)
deallocate(alfa)
deallocate(vector)
deallocate(all_det)

return
end subroutine vector_rep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine overlap(nl,str1,totstr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
implicit none

integer::i,i1,i2,i3,i4,i5,i6,i7,j,k,l,col1(2,1000),nl,str1(15000,20),totstr,nio
real::orb_ovlp_mat1_norm(10,10)
Double Precision::D1(1000),aldet_new
double precision::ovlp_mat(5000,100)

nio=niao
if(ovlp_int.eq.0)then
  do i=1,nactorb
    do j=1,nactorb
      if(i.eq.j)then
      orb_ovlp_mat1(i,j)=4.0
      else
        if(atm_nb_sym(i).eq.atm_nb_sym(j))then 
          if(atm_nb_orbs(i).eq.atm_nb_orbs(j))then
          orb_ovlp_mat1(i,j)=3.0
          else
          orb_ovlp_mat1(i,j)=2.0/dist_rel_mat(atm_nb_orbs(i),atm_nb_orbs(j))
          endif
        else
        orb_ovlp_mat1(i,j)=0.0
        endif
      endif
    enddo
  enddo
endif

do i1=1,5000
  do i2=1,100
    ovlp_mat(i1,i2)=0.0
  enddo
enddo


do i1=1,totstr
  do i2=i1,totstr
    do i3=1,ndet
      do i4=1,ndet
      aldet_new=1.0
        do i5=1,nalpha+nl
        col1(1,i5)=detmnt(str_det_sec(i1,i3),i5)-nio
        enddo
        do i5=1,nalpha+nl
        col1(2,i5)=detmnt(str_det_sec(i2,i4),i5)-nio
        enddo
        call MatLDR1('orb',col1,nalpha+nl,D1)
        do i6=1,nalpha+nl
        aldet_new=aldet_new*D1(i6)
        enddo
        do i6=1,20
        col1(1,i6)=0
        col1(2,i6)=0
        enddo
        i6=0
        do i5=nalpha+nl+1,nalpha+nl+nbeta+nl
        i6=i6+1
        col1(1,i6)=detmnt(str_det_sec(i1,i3),i5)-nio
        enddo
        i6=0
        do i5=nalpha+nl+1,nalpha+nl+nbeta+nl
        i6=i6+1
        col1(2,i6)=detmnt(str_det_sec(i2,i4),i5)-nio
        enddo
        call MatLDR1('orb',col1,nbeta+nl,D1)
        do i6=1,nbeta+nl
        aldet_new=aldet_new*D1(i6)
        enddo
        ovlp_mat(i1,i2)=ovlp_mat(i1,i2)+aldet_new*det_sign(str_det_sec(i1,i3))*det_sign(str_det_sec(i2,i4))
      enddo
    enddo
  enddo
enddo

do i1=1,totstr
  do i2=i1+1,totstr
  ovlp_mat(i2,i1)=ovlp_mat(i1,i2)
  enddo
enddo

do i1=1,totstr
  do i2=1,totstr
  ovlp_mat_norm(i1,i2)=ovlp_mat(i1,i2)/sqrt(ovlp_mat(i1,i1)*ovlp_mat(i2,i2))
  enddo
enddo

321 format(I5,a,50I3)
320 format(I5,a,50F11.3)

return
end subroutine overlap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine matrx_norm(ovlp_mat,totstr,ovlp_mat_norm)
implicit none

integer::i1,i2,totstr
double precision::ovlp_mat(5000,100),ovlp_mat_norm(5000,100)

do i1=1,totstr
do i2=1,totstr
ovlp_mat_norm(i1,i2)=ovlp_mat(i1,i2)/sqrt(ovlp_mat(i1,i1)*ovlp_mat(i2,i2))
enddo
enddo


end subroutine matrx_norm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine intra_bond_factor(nl,str,tonstruc,str_quality_1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::nl,i2,m8,l1,l2,l3,k1,k2,m13,m14,str(15000,20),str_quality_1(15000),tonstruc

print*,'enter intra_bond_factor'
do m8=1,15000
str_quality_1(m8)=1
enddo

l3=1
do m8=1,tonstruc
!print*,'ssttr',(str(m8,i2),i2=1,nae)
l2=1
do k2=1+nl*2,nae-nlast,2
do m13=1,atom
l1=0
do m14=1,atn(active_atoms(m13))
if(atn(active_atoms(m13)).eq.1)goto 505
do k1=k2,k2+1
if(str(m8,k1).eq.atoset(active_atoms(m13),m14))then
l1=l1+1
endif
enddo
505 enddo
if(l1.eq.2) then
l2=l2+l3
goto 506
endif
enddo
506 enddo
str_quality_1(m8)=l2
!print*,'intra_bond_factor',str_quality_1(m8)
enddo


print*,'exit intra_bond_factor'
return
end subroutine intra_bond_factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine symm_break_factor(nl,str,tonstruc,str_quality_2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::nl,i2,m8,l1,l2,l3,k1,k2,m13,m14,str(15000,20),str_quality_2(15000),tonstruc
integer::atsymset(20,20),nsym,syn(50),at_sym(50)

common/ats/atsymset,nsym,syn,at_sym

print*,'enter symm_break_factor'
do m8=1,15000
str_quality_2(m8)=1
enddo
do m8=1,tonstruc
!print*,'symm_break:sttr',(str(m8,i2),i2=1,nae)
l2=1
if(nsym.eq.1)goto 509
if(nsym.ne.1)l2=1+(nae-nlast-nl*2)/2
do k2=1+nl*2,nae-nlast,2
!print*,'k2',k2
do m13=1,nsym
l1=0
do m14=1,syn(m13)
if(syn(m13).eq.1)goto 507
do k1=k2,k2+1
if(str(m8,k1).eq.atsymset(m13,m14))then
l1=l1+1
endif
enddo
507 enddo
if(l1.eq.2) then
l2=l2-1 
goto 508
endif
enddo
508 enddo
str_quality_2(m8)=l2
509 enddo

!print*,'tonstruc',tonstruc
do m8=1,tonstruc
!!quality_fac(m8)=str_quality_1(m8)*str_quality_2(m8)-(str_quality_2(m8)-1)*(str_quality_1(m8)-1)
print*,'symm_break:sttr',(str(m8,i2),i2=1,nae)
print*,'qqqqqq',str_quality_2(m8)
enddo

print*,'exit symm_break_factor'
return
end subroutine symm_break_factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!subroutine w_s_x(nl,strn,str3,ncqs,q_fac2)
!!!
!!!!!!! this subroutine select the independent structures and write them in the .xmi file 
!!!!!!! for the symmetric set of structures of a symmetric molecule (as opted by user) ....
!!!!!!! subroutine started to select the structures of highest quality biggest
!!!!!!! symmetric set and select the structures as a set of symmetric set if
!!!!!!! they are independent with the previous sets or reject the full set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!use commondat
!!implicit none
!!common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
!!sigsym,tnqs_sig
!!common/str/str5,nstr7
!!
!!integer::nl,strn,ncqs,tostr,initstr,i,i1,i2,i3,i4,i5,i6,i7,i8,i9,m119,m18,m19,m20,m21,m23,m24,count&
!!,qul(100),nqul,j,jj,jjj,fg,flg,ii5,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,x,y,strs
!!integer::i5up16,i16i7,m16m21,i5up15,i15i7,m15m21,i5up14,i14i7,m14m21,i5up13,i13i7,m13m21,i5up12,i12i7,m12m21,&
!!i5up11,i11i7,m11m21,i5up10,i10i7,m10m21,i5up9,i9i7,m9m21,i5up8,i8i7,m8m21,i5up7,i7i7,m7m21,i5up6,i6i7,m6m21,&
!!i5up5,i5i7,m5m21,i5up4,i4i7,m4m21,i5up3,i3i7,m3m21,i5up2,i2i7,m2m21,i5up1,i1i7,m1m21,i5up17,i17i7,m17m21
!!integer::str3(15000,20),q_fac2(15000),finalvec(15000),strset(1000),col(1000),sigsym(15000),tnqs_sig,&
!!ffvec2(15000,1000),bondq(15000),bondq4(15000),nqset(15000),str5(2000,20),nstr7,&
!!tndet,totstr,Ifail,indpnt,strno(1000),str_quality_1(15000),str_quality_2(15000),ttqlty0,ttqlty&
!!,tqlty,bqlty,sqlty,hqlty,tnqs,nssym,qulsym(15000),symq(15000),set_number,ttqlty1,det_inv,ttqlty2,ttqlty3
!!integer::rumer(15000),rumer_rad(15000),quality_fac(15000),rumset,u1,max_set,mns
!!!integer,dimension(:),allocatable::qq1,qq2,qq
!!integer::qq1(5000),qq2(5000),qq(5000),str2(2000,20),Rid,set_num(100)
!!real*8::ovlp
!!Double Precision::D(1000)
!!character(10)::dd,a
!!character(len=100)::outfile
!!
!!
!!print*,'enter write_symm_xmi'
!!!x=100000
!!!y=100
!!!
!!!allocate(str2(x,y))
!!!allocate(qq(x))
!!!allocate(qq1(x))
!!!allocate(qq2(x))
!!!str2(x,y)=0
!!!qq(x)=0
!!!qq1(x)=0
!!!qq2(x)=0
!!mns=0
!!u1=1
!!max_set=75000
!!
!!if(nfset.eq.3.or.nfset.eq.5)then
!!rumset=0
!!call rumer_structures(nl,str3,ncqs,rumer,rumer_rad)
!!call write_rumer_xmi(nl,str3,ncqs,rumer,rumer_rad,quality_fac)
!!endif
!!
!!!print *,'write_symm_xmi'
!!set_number=0
!!bqlty=0
!!tqlty=0
!!sqlty=0
!!hqlty=0
!!indpnt=2
!!!ovlpval=1.0
!!if(noq0.gt.strn)then
!!ttqlty0=noq0
!!else
!!ttqlty0=strn+noq0
!!endif
!!ttqlty1=strn+noq1
!!ttqlty2=strn+noq2
!!ttqlty3=strn+noq3
!!
!!write(*,*)'ttqlty1,ttqlty2,ttqlty3',ttqlty1,ttqlty2,ttqlty3,ttqlty0,noq1,noq2,noq3,noq0,nfset
!!!!!!!!! strn = number of permissible structures !!!!
!!!print*,'sou1',strn,ndet,nl,ncqs
!!do i=1,ncqs
!!write(*,231),i,(str3(i,j),j=1,nae),q_fac2(i),str_quality_1(i),str_quality_2(i),bondq(i)
!!enddo
!!!if(indpnt.eq.1)then
!!!call vector_rep(nl,str3,ncqs,ffvec2)
!!!endif
!!!print*,'sou1',strn,ndet,nl,ncqs
!!!do i1=1,ncqs
!!!print*,(ffvec2(i1,i),i=1,ndet)
!!!enddo
!!!call quality_factor(nl,str3,ncqs,q_fac2,q_1,q_2)
!!jj=1
!!do m19=1,ncqs
!!!print*,'q_fac2',q_fac2(m19)
!!if(m19.eq.1)qul(1)=q_fac2(1)
!!j=jj
!!do i=1,j
!!if(qul(i).eq.q_fac2(m19))goto 373
!!enddo
!!jj=jj+1
!!qul(i)=q_fac2(m19)
!!!print*,qul(i)
!!373 enddo
!!nqul=jj
!!!print*,'ncqs',ncqs
!!
!!do i=1,nqul
!!jjj=0
!!jj=0
!!do m19=1,ncqs
!!!print*,qul(i),q_fac2(m19)
!!if(qul(i).eq.q_fac2(m19))then
!!jjj=jjj+1
!!jj=m19
!!endif
!!enddo
!!nqset(i)=jjj
!!strset(i)=jj
!!print*,'nqset(i)',i,nqset(i),qul(i),strset(i)
!!!write(*,231),(str3(i,j),j=1,nae)
!!enddo
!!
!!flg=0
!!totstr=0
!!i7=0
!!m21=0
!!i5=ndet
!!do i8=1,1000
!!finalvec(i8)=0
!!enddo
!!do m19=1,10000
!!do m18=1,15
!!str2(m19,m18)=0
!!enddo
!!enddo
!!
!!
!!i5up1=i5
!!i1i7=i7
!!m1m21=m21
!!
!!!!!!! LOOP STARTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!!!!! loop 1 >>>
!!do m1=1,nqul
!!print*,'loop1'
!!
!!do i=i5up1+1,i5
!!finalvec(i)=0
!!enddo
!!
!!do m19=i1i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!
!!i7=i1i7
!!totstr=i1i7
!!m21=m1m21
!!
!!jj=0
!!
!!flg=0
!!
!!if(nqset(m1).gt.strn)goto 701
!!
!!if(m1.eq.1)then
!!do i=1,strset(m1)
!!totstr=totstr+1
!!strno(totstr)=i
!!enddo
!!else
!!do i=strset(m1-1)+1,strset(m1)
!!totstr=totstr+1
!!strno(totstr)=i
!!enddo
!!endif
!!if(totstr.gt.strn)goto 701
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!if(Ifail.eq.1)goto 701
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!call MatLDR('str',strno,totstr,D)
!!
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!if(1.0-ovlp.gt.ovval)goto 701
!!endif
!!
!!m21=m21+nqset(m1)
!!if(m21.gt.strn) goto 701
!!
!!if(m1.eq.1)strs=0
!!if(m1.ne.1)strs=strset(m1-1)
!!
!!do m19=strs+1,strset(m1)
!!if(q_fac2(m19).ne.qul(m1))goto 601
!!i7=i7+1
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=1
!!601 enddo
!!if(i7.eq.strn) then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!enddo
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!
!!goto 701
!!endif
!!
!!i5up2=i5
!!i2i7=i7
!!m2m21=m21
!!
!!!!!!! loop 2 >>>
!!
!!do m2=m1+1,nqul
!!print*,'loop2'
!!
!!do i=i5up2+1,i5
!!finalvec(i)=0
!!enddo
!!
!!do m19=i2i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!
!!i7=i2i7
!!totstr=i2i7
!!m21=m2m21
!!
!!jj=0
!!
!!flg=0
!!
!!if(nqset(m2).gt.strn)goto 702
!!
!!do i=strset(m2-1)+1,strset(m2)
!!totstr=totstr+1
!!strno(totstr)=i
!!enddo
!!if(totstr.gt.strn)goto 702
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!
!!if(Ifail.eq.1) goto 702
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!call MatLDR('str',strno,totstr,D)
!!
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!if(1.0-ovlp.gt.ovval)goto 702
!!endif
!!
!!
!!m21=m21+nqset(m2)
!!if(m21.gt.strn) goto 702
!!
!!
!!do m19=strset(m2-1)+1,strset(m2)
!!if(q_fac2(m19).ne.qul(m2))goto 602
!!i7=i7+1
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=2
!!602 enddo
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!enddo
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!goto 702
!!endif
!!
!!
!!i5up3=i5
!!i3i7=i7
!!m3m21=m21
!!
!!!!!!! loop 3 >>>
!!
!!do m3=m2+1,nqul
!!print*,'loop3'
!!
!!do i=i5up3+1,i5
!!finalvec(i)=0
!!enddo
!!
!!do m19=i3i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!
!!i7=i3i7
!!totstr=i3i7
!!m21=m3m21
!!
!!jj=0
!!
!!flg=0
!!
!!!do m19=strset(m3-1)+1,strset(m3)
!!
!!if(nqset(m3).gt.strn)goto 703
!!
!!do i=strset(m3-1)+1,strset(m3)
!!totstr=totstr+1
!!strno(totstr)=i
!!print*,'strno',strno(totstr)
!!enddo
!!if(totstr.gt.strn)goto 703
!!do i=1,totstr
!!print*,'strno_all',strno(i)
!!enddo
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!
!!print*,'Ifail',Ifail,i7,q_fac2(m19)
!!if(Ifail.eq.1) goto 703
!!!do i=strset(m3-1)+1,strset(m3)
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!enddo
!!!goto 703
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 703
!!endif
!!
!!m21=m21+nqset(m3)
!!if(m21.gt.strn) goto 703
!!!print*,'sourav_xmi_sym'
!!
!!!303 enddo
!!
!!
!!do m19=strset(m3-1)+1,strset(m3)
!!if(q_fac2(m19).ne.qul(m3))goto 603
!!i7=i7+1
!!print*,'new struc',i7,m19,q_fac2(m19),str_quality_1(m19),str_quality_2(m19),bondq(m19)
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=3
!!603 enddo
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!print*,'ttqlty,ttqlty0',ttqlty,ttqlty0
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!!endif
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(ovopt.eq.vpt) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!goto 703
!!endif
!!
!!i5up4=i5
!!i4i7=i7
!!m4m21=m21
!!
!!!!! loop 4 >>>
!!do m4=m3+1,nqul
!!print*,'loop4'
!!!print*,'m4m4m4',m4
!!do i=i5up4+1,i5
!!finalvec(i)=0
!!enddo
!!do m19=i4i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!i7=i4i7
!!totstr=i4i7
!!m21=m4m21
!!
!!jj=0
!!
!!flg=0
!!!do m19=strset(m4-1)+1,strset(m4)
!!
!!if(nqset(m4).gt.strn)goto 704
!!
!!do i=strset(m4-1)+1,strset(m4)
!!totstr=totstr+1
!!strno(totstr)=i
!!print*,'strno',strno(totstr)
!!enddo
!!
!!if(totstr.gt.strn)goto 704
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!!print*,'ifail4',ifail
!!print*,'Ifail',Ifail,i7,q_fac2(m19)
!!if(Ifail.eq.1) goto 704
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!goto 304
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 704
!!endif
!!
!!
!!m21=m21+nqset(m4)
!!if(m21.gt.strn) goto 704
!!!print*,'sourav_xmi_sym'
!!
!!!304 enddo
!!
!!
!!do m19=strset(m4-1)+1,strset(m4)
!!if(q_fac2(m19).ne.qul(m4))goto 604
!!i7=i7+1
!!print*,'new struc',i7,m19,q_fac2(m19),str_quality_1(m19),str_quality_2(m19),bondq(m19)
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=4
!!604 enddo
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!print*,'ttqlty,ttqlty0',ttqlty,ttqlty0
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!!endif
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(ovopt.eq.vpt) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!goto 704
!!endif
!!
!!i5up5=i5
!!i5i7=i7
!!m5m21=m21
!!
!!!!!!! loop 5 >>>
!!do m5=m4+1,nqul
!!print*,'loop5'
!!!print*,'m5m5m5',m5
!!do i=i5up5+1,i5
!!finalvec(i)=0
!!enddo
!!do m19=i5i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!i7=i5i7
!!totstr=i5i7
!!!print*,'totstr5',totstr
!!m21=m5m21
!!
!!jj=0
!!
!!flg=0
!!!do m19=strset(m5-1)+1,strset(m5)
!!
!!if(nqset(m5).gt.strn)goto 705
!!
!!do i=strset(m5-1)+1,strset(m5)
!!totstr=totstr+1
!!strno(totstr)=i
!!print*,'strno',strno(totstr)
!!enddo
!!
!!if(totstr.gt.strn)goto 705
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!print*,'Ifail',Ifail,i7,q_fac2(m19)
!!if(Ifail.eq.1) goto 705
!!
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!goto 305
!!!endif
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 705
!!endif
!!
!!
!!m21=m21+nqset(m5)
!!if(m21.gt.strn) goto 705
!!!print*,'sourav_xmi_sym'
!!
!!!305 enddo
!!
!!!print*,'jj',jj,m5,nqul
!!!print*,'******jj',jj
!!
!!do m19=strset(m5-1)+1,strset(m5)
!!if(q_fac2(m19).ne.qul(m5))goto 605
!!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!!i7=i7+1
!!print*,'new struc',i7,m19,q_fac2(m19),str_quality_1(m19),str_quality_2(m19),bondq(m19)
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=5
!!605 enddo
!!!print*,'i7 5',i7,m5
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!print*,'ttqlty,ttqlty0',ttqlty,ttqlty0
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!!endif
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!goto 705
!!endif
!!
!!i5up6=i5
!!i6i7=i7
!!m6m21=m21
!!
!!do m6=m5+1,nqul
!!
!!print*,'loop6'
!!!print*,'m6m6m6',m6
!!do i=i5up6+1,i5
!!finalvec(i)=0
!!enddo
!!do m19=i6i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!i7=i6i7
!!totstr=i6i7
!!!print*,'totstr6',totstr
!!m21=m6m21
!!
!!jj=0
!!
!!flg=0
!!!do m19=strset(m6-1)+1,strset(m6)
!!
!!if(nqset(m6).gt.strn)goto 706
!!
!!do i=strset(m6-1)+1,strset(m6)
!!totstr=totstr+1
!!strno(totstr)=i
!!print*,'strno',strno(totstr)
!!enddo
!!
!!if(totstr.gt.strn)goto 706
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!!print*,'ifail6',ifail
!!print*,'Ifail',Ifail,i7,q_fac2(m19),str_quality_1(m19),str_quality_2(m19),bondq(m19)
!!if(Ifail.eq.1)goto 706
!!!!do i8=1,nae
!!!!str4(totstr,i8)=0
!!!!enddo
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!goto 306
!!!endif
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 706
!!endif
!!
!!
!!m21=m21+nqset(m6)
!!if(m21.gt.strn) goto 706
!!!print*,'sourav_xmi_sym'
!!
!!!306 enddo
!!
!!!print*,'jj',jj,m6,nqul
!!!print*,'******jj',jj
!!
!!do m19=strset(m6-1)+1,strset(m6)
!!if(q_fac2(m19).ne.qul(m6))goto 606
!!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!!i7=i7+1
!!print*,'new struc',i7,m19,q_fac2(m19),str_quality_1(m19),str_quality_2(m19),bondq(m19)
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=6
!!606 enddo
!!!print*,'i7  6',i7,m6
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!print*,'ttqlty,ttqlty0',ttqlty,ttqlty0
!!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0,tqlty,ttqlty1,bqlty,ttqlty2,sqlty,ttqlty3
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!!endif
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!goto 706
!!endif
!!
!!i5up7=i5
!!i7i7=i7
!!m7m21=m21
!!
!!do m7=m6+1,nqul
!!print*,'loop7'
!!!print*,'m7m7m7',m7
!!do i=i5up7+1,i5
!!finalvec(i)=0
!!enddo
!!do m19=i7i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!i7=i7i7
!!totstr=i7i7
!!!print*,'totstr7',totstr
!!m21=m7m21
!!
!!jj=0
!!
!!flg=0
!!!do m19=strset(m7-1)+1,strset(m7)
!!
!!if(nqset(m7).gt.strn)goto 707
!!
!!do i=strset(m7-1)+1,strset(m7)
!!totstr=totstr+1
!!strno(totstr)=i
!!print*,'strno',strno(totstr)
!!enddo
!!if(totstr.gt.strn)goto 707
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!print*,'Ifail',Ifail,i7,q_fac2(m19)
!!if(Ifail.eq.1)goto 707
!!!!do i8=1,nae
!!!!str4(totstr,i8)=0
!!!!enddo
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!goto 307
!!!endif
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 707
!!endif
!!
!!
!!m21=m21+nqset(m7)
!!if(m21.gt.strn) goto 707
!!!print*,'sourav_xmi_sym'
!!
!!!307 enddo
!!
!!!print*,'jj',jj,m7,nqul
!!!407 if(jj.eq.nqset(m7))then
!!!print*,'******jj',jj
!!
!!do m19=strset(m7-1)+1,strset(m7)
!!if(q_fac2(m19).ne.qul(m7))goto 607
!!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!!i7=i7+1
!!print*,'new struc',i7,m19,q_fac2(m19),str_quality_1(m19),str_quality_2(m19),bondq(m19)
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=7
!!607 enddo
!!print*,'i7 7',i7,strn
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!print*,'ttqlty,ttqlty0',ttqlty,ttqlty0
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!!endif
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!goto 707
!!endif
!!
!!i5up8=i5
!!i8i7=i7
!!m8m21=m21
!!
!!do m8=m7+1,nqul
!!print*,'loop8'
!!!print*,'m8m8m8',m8
!!do i=i5up8+1,i5
!!finalvec(i)=0
!!enddo
!!do m19=i8i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!i7=i8i7
!!totstr=i8i7
!!!print*,'totstr8',totstr
!!m21=m8m21
!!
!!jj=0
!!
!!flg=0
!!!do m19=strset(m8-1)+1,strset(m8)
!!
!!if(nqset(m8).gt.strn)goto 708
!!
!!do i=strset(m8-1)+1,strset(m8)
!!totstr=totstr+1
!!strno(totstr)=i
!!print*,'strno',strno(totstr)
!!enddo
!!
!!if(totstr.gt.strn)goto 708
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!print*,'Ifail',Ifail,i7,q_fac2(m19)
!!if(Ifail.eq.1)goto 708
!!
!!!!do i8=1,nae
!!!!str4(totstr,i8)=0
!!!!enddo
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!goto 308
!!!endif
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 708
!!endif
!!
!!
!!m21=m21+nqset(m8)
!!if(m21.gt.strn) goto 708
!!!print*,'sourav_xmi_sym'
!!
!!!308 enddo
!!
!!!print*,'jj',jj,m7,nqul
!!!408 if(jj.eq.nqset(m8))then
!!!print*,'******jj',jj
!!
!!do m19=strset(m8-1)+1,strset(m8)
!!if(q_fac2(m19).ne.qul(m8))goto 608
!!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!!i7=i7+1
!!print*,'new struc',i7,m19,q_fac2(m19)
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=8
!!608 enddo
!!!print*,'i7 8',i7,m8
!!print*,'i7 7',i7,strn
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!!endif
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!goto 708
!!endif
!!
!!i5up9=i5
!!i9i7=i7
!!m9m21=m21
!!
!!do m9=m8+1,nqul
!!print*,'loop9'
!!!print*,'m9m9m9',m9
!!do i=i5up9+1,i5
!!finalvec(i)=0
!!enddo
!!do m19=i9i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!i7=i9i7
!!totstr=i9i7
!!!print*,'totstr9',totstr
!!m21=m9m21
!!
!!jj=0
!!
!!flg=0
!!!do m19=strset(m9-1)+1,strset(m9)
!!
!!if(nqset(m9).gt.strn)goto 709
!!
!!do i=strset(m9-1)+1,strset(m9)
!!totstr=totstr+1
!!strno(totstr)=i
!!print*,'strno',strno(totstr)
!!enddo
!!
!!if(totstr.gt.strn)goto 709
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!!print*,'ifail9',ifail
!!print*,'Ifail',Ifail,i7,q_fac2(m19)
!!if(Ifail.eq.1)goto 709
!!!!do i8=1,nae
!!!!str4(totstr,i8)=0
!!!!enddo
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!goto 309
!!!endif
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 709
!!endif
!!
!!
!!m21=m21+nqset(m9)
!!if(m21.gt.strn) goto 709
!!!print*,'sourav_xmi_sym'
!!
!!!309 enddo
!!
!!!print*,'jj',jj,m9,nqul
!!!409 if(jj.eq.nqset(m9))then
!!!print*,'******jj',jj
!!
!!do m19=strset(m9-1)+1,strset(m9)
!!if(q_fac2(m19).ne.qul(m9))goto 609
!!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!!i7=i7+1
!!print*,'new struc',i7,m19,q_fac2(m19)
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=9
!!609 enddo
!!!print*,'i7 9',i7,m9
!!print*,'i7 7',i7,strn
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!!endif
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!goto 709
!!endif
!!
!!i5up10=i5
!!i10i7=i7
!!m10m21=m21
!!
!!do m10=m9+1,nqul
!!print*,'loop10'
!!!print*,'m10m10m10',m10
!!do i=i5up10+1,i5
!!finalvec(i)=0
!!enddo
!!do m19=i10i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!i7=i10i7
!!totstr=i10i7
!!!print*,'totstr10',totstr
!!m21=m10m21
!!
!!jj=0
!!
!!flg=0
!!!do m19=strset(m10-1)+1,strset(m10)
!!
!!if(nqset(m10).gt.strn)goto 710
!!
!!do i=strset(m10-1)+1,strset(m10)
!!totstr=totstr+1
!!strno(totstr)=i
!!print*,'strno',strno(totstr)
!!enddo
!!
!!if(totstr.gt.strn)goto 710
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!print*,'Ifail',Ifail,i7,q_fac2(m19)
!!if(Ifail.eq.1)goto 710
!!!!do i8=1,nae
!!!!str4(totstr,i8)=0
!!!!enddo
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!goto 310
!!!endif
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 710
!!endif
!!
!!
!!m21=m21+nqset(m10)
!!if(m21.gt.strn) goto 710
!!!print*,'sourav_xmi_sym'
!!
!!!310 enddo
!!
!!!print*,'jj',jj,m10,nqul
!!!410 if(jj.eq.nqset(m10))then
!!!print*,'******jj',jj
!!
!!do m19=strset(m10-1)+1,strset(m10)
!!if(q_fac2(m19).ne.qul(m10))goto 610
!!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!!i7=i7+1
!!print*,'new struc',i7,m19
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=10
!!610 enddo
!!!print*,'i7 10',i7,m10
!!print*,'i7 7',i7,strn
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!!endif
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!goto 710
!!endif
!!
!!i5up11=i5
!!i11i7=i7
!!m11m21=m21
!!
!!do m11=m10+1,nqul
!!print*,'loop11'
!!!print*,'m11m11m11',m11
!!do i=i5up11+1,i5
!!finalvec(i)=0
!!enddo
!!do m19=i11i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!i7=i11i7
!!totstr=i11i7
!!!print*,'totstr11',totstr
!!m21=m11m21
!!
!!jj=0
!!
!!flg=0
!!!do m19=strset(m11-1)+1,strset(m11)
!!
!!if(nqset(m11).gt.strn)goto 711
!!
!!do i=strset(m11-1)+1,strset(m11)
!!totstr=totstr+1
!!strno(totstr)=i
!!enddo
!!
!!if(totstr.gt.strn)goto 711
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!print*,'Ifail',Ifail,i7
!!if(Ifail.eq.1)goto 711
!!!!do i8=1,nae
!!!!str4(totstr,i8)=0
!!!!enddo
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!goto 311
!!!endif
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 711
!!endif
!!
!!
!!m21=m21+nqset(m11)
!!if(m21.gt.strn) goto 711
!!!print*,'sourav_xmi_sym'
!!
!!!311 enddo
!!
!!!print*,'jj',jj,m11,nqul
!!!411 if(jj.eq.nqset(m11))then
!!!print*,'******jj',jj
!!
!!do m19=strset(m11-1)+1,strset(m11)
!!if(q_fac2(m19).ne.qul(m11))goto 611
!!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!!i7=i7+1
!!print*,'new struc',i7,m19
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=11
!!611 enddo
!!!print*,'i7 11',i7,m11
!!print*,'i7 7',i7,strn
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!!endif
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!goto 711
!!endif
!!
!!i5up12=i5
!!i12i7=i7
!!m12m21=m21
!!
!!do m12=m11+1,nqul
!!print*,'loop12'
!!!print*,'m12m12m12',m12
!!do i=i5up12+1,i5
!!finalvec(i)=0
!!enddo
!!do m19=i12i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!i7=i12i7
!!totstr=i12i7
!!!print*,'totstr12',totstr
!!m21=m12m21
!!
!!jj=0
!!
!!flg=0
!!!do m19=strset(m12-1)+1,strset(m12)
!!
!!if(nqset(m12).gt.strn)goto 712
!!
!!do i=strset(m12-1)+1,strset(m12)
!!totstr=totstr+1
!!strno(totstr)=i
!!enddo
!!
!!if(totstr.gt.strn)goto 712
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!print*,'Ifail',Ifail,i7
!!if(Ifail.eq.1)goto 712
!!!!do i8=1,nae
!!!!str4(totstr,i8)=0
!!!!enddo
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!goto 312
!!!endif
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 712
!!endif
!!
!!
!!m21=m21+nqset(m12)
!!if(m21.gt.strn) goto 712
!!!print*,'sourav_xmi_sym'
!!
!!
!!!312 enddo
!!
!!!print*,'jj',jj,m12,nqul
!!!412 if(jj.eq.nqset(m12))then
!!!print*,'******jj',jj
!!
!!do m19=strset(m12-1)+1,strset(m12)
!!if(q_fac2(m19).ne.qul(m12))goto 612
!!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!!i7=i7+1
!!print*,'new struc',i7,m19
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=12
!!612 enddo
!!!print*,'i7 12',i7,m12
!!print*,'i7 7',i7,strn
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!!endif
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!goto 712
!!endif
!!
!!i5up13=i5
!!i13i7=i7
!!m13m21=m21
!!
!!do m13=m12+1,nqul
!!print*,'loop13'
!!
!!!print*,'m13m13m13',m13
!!do i=i5up13+1,i5
!!finalvec(i)=0
!!enddo
!!do m19=i13i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!i7=i13i7
!!totstr=i13i7
!!!print*,'totstr13',totstr
!!m21=m13m21
!!
!!jj=0
!!
!!flg=0
!!!do m19=strset(m13-1)+1,strset(m13)
!!
!!if(nqset(m13).gt.strn)goto 713
!!
!!do i=strset(m13-1)+1,strset(m13)
!!totstr=totstr+1
!!strno(totstr)=i
!!enddo
!!
!!if(totstr.gt.strn)goto 713
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!!print*,'ifail13',ifail
!!print*,'Ifail',Ifail,i7
!!if(Ifail.eq.1)goto 713
!!!!do i8=1,nae
!!!!str4(totstr,i8)=0
!!!!enddo
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!goto 313
!!!endif
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 713
!!endif
!!
!!m21=m21+nqset(m13)
!!if(m21.gt.strn) goto 713
!!!print*,'sourav_xmi_sym'
!!
!!!313 enddo
!!
!!!print*,'jj',jj,m13,nqul
!!!413 if(jj.eq.nqset(m13))then
!!!print*,'******jj',jj
!!
!!do m19=strset(m13-1)+1,strset(m13)
!!if(q_fac2(m19).ne.qul(m13))goto 613
!!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!!i7=i7+1
!!print*,'new struc',i7,m19
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=13
!!613 enddo
!!!print*,'i7 13',i7,m13
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!goto 713
!!endif
!!
!!i5up14=i5
!!i14i7=i7
!!m14m21=m21
!!
!!do m14=m13+1,nqul
!!print*,'loop14'
!!!print*,'m14m14m14',m14
!!do i=i5up14+1,i5
!!finalvec(i)=0
!!enddo
!!do m19=i14i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!i7=i14i7
!!totstr=i14i7
!!!print*,'totstr14',totstr
!!m21=m14m21
!!
!!jj=0
!!
!!flg=0
!!!do m19=strset(m14-1)+1,strset(m14)
!!
!!if(nqset(m14).gt.strn)goto 714
!!
!!do i=strset(m14-1)+1,strset(m14)
!!totstr=totstr+1
!!strno(totstr)=i
!!enddo
!!
!!if(totstr.gt.strn)goto 714
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!print*,'Ifail',Ifail,i7
!!if(Ifail.eq.1)goto 714
!!!!do i8=1,nae
!!!!str4(totstr,i8)=0
!!!!enddo
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!goto 314
!!!endif
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 714
!!endif
!!
!!
!!m21=m21+nqset(m14)
!!if(m21.gt.strn) goto 714
!!!print*,'sourav_xmi_sym'
!!
!!!314 enddo
!!
!!!print*,'jj',jj,m14,nqul
!!!414 if(jj.eq.nqset(m14))then
!!!print*,'******jj',jj
!!
!!do m19=strset(m14-1)+1,strset(m14)
!!if(q_fac2(m19).ne.qul(m14))goto 614
!!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!!i7=i7+1
!!print*,'new struc',i7,m19
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=14
!!614 enddo
!!!print*,'i7 14',i7,m14
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!!endif
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!goto 714
!!endif
!!
!!i5up15=i5
!!i15i7=i7
!!m15m21=m21
!!
!!do m15=m14+1,nqul
!!print*,'loop15'
!!!print*,'m15m15m15',m15
!!do i=i5up15+1,i5
!!finalvec(i)=0
!!enddo
!!do m19=i15i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!i7=i15i7
!!totstr=i15i7
!!!print*,'totstr15',totstr
!!m21=m15m21
!!
!!jj=0
!!
!!flg=0
!!!do m19=strset(m15-1)+1,strset(m15)
!!
!!if(nqset(m15).gt.strn)goto 715
!!
!!do i=strset(m15-1)+1,strset(m15)
!!totstr=totstr+1
!!strno(totstr)=i
!!enddo
!!
!!if(totstr.gt.strn)goto 715
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!print*,'Ifail',Ifail,i7
!!if(Ifail.eq.1)goto 715
!!!!do i8=1,nae
!!!!str4(totstr,i8)=0
!!!!enddo
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!goto 315
!!!endif
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 715
!!endif
!!
!!
!!m21=m21+nqset(m15)
!!if(m21.gt.strn) goto 715
!!!print*,'sourav_xmi_sym'
!!
!!!315 enddo
!!
!!!print*,'jj',jj,m15,nqul
!!!415 if(jj.eq.nqset(m15))then
!!!print*,'******jj',jj
!!
!!do m19=strset(m15-1)+1,strset(m15)
!!if(q_fac2(m19).ne.qul(m15))goto 615
!!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!!i7=i7+1
!!print*,'new struc',i7,m19
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=15
!!615 enddo
!!!print*,'i7 15',i7,m15
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!!endif
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!goto 715
!!endif
!!
!!i5up16=i5
!!i16i7=i7
!!m16m21=m21
!!
!!do m16=m15+1,nqul
!!print*,'loop16'
!!!print*,'m16m16m16',m16
!!do i=i5up16+1,i5
!!finalvec(i)=0
!!enddo
!!do m19=i16i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!i7=i16i7
!!totstr=i16i7
!!!print*,'totstr16',totstr
!!m21=m16m21
!!
!!jj=0
!!
!!flg=0
!!!do m19=strset(m16-1)+1,strset(m16)
!!
!!if(nqset(m16).gt.strn)goto 716
!!
!!do i=strset(m16-1)+1,strset(m16)
!!totstr=totstr+1
!!strno(totstr)=i
!!enddo
!!if(totstr.gt.strn)goto 716
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!print*,'Ifail',Ifail,i7
!!if(Ifail.eq.1)goto 716
!!!!do i8=1,nae
!!!!str4(totstr,i8)=0
!!!!enddo
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!goto 316
!!!endif
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 716
!!endif
!!
!!
!!m21=m21+nqset(m16)
!!if(m21.gt.strn) goto 716
!!!print*,'sourav_xmi_sym'
!!
!!!316 enddo
!!
!!!print*,'jj',jj,m16,nqul
!!!416 if(jj.eq.nqset(m16))then
!!!print*,'******jj',jj
!!
!!do m19=strset(m16-1)+1,strset(m16)
!!if(q_fac2(m19).ne.qul(m16))goto 616
!!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!!i7=i7+1
!!print*,'new struc',i7,m19
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!count=16
!!616 enddo
!!!print*,'i7 16',i7,m16
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!!endif
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!goto 716
!!endif
!!
!!
!!i5up17=i5
!!i17i7=i7
!!m17m21=m21
!!
!!do m17=m16+1,nqul
!!print*,'loop17'
!!!print*,'m17m17m17',m17
!!do i=i5up17+1,i5
!!finalvec(i)=0
!!enddo
!!do m19=i17i7+1,i7
!!strno(m19)=0
!!do i8=1,15
!!str2(m19,i8)=0
!!enddo
!!enddo
!!i7=i17i7
!!totstr=i17i7
!!!print*,'totstr17',totstr
!!m21=m17m21
!!
!!jj=0
!!
!!flg=0
!!!do m19=strset(m17-1)+1,strset(m17)
!!
!!if(nqset(m17).gt.strn)goto 717
!!
!!do i=strset(m17-1)+1,strset(m17)
!!totstr=totstr+1
!!strno(totstr)=i
!!enddo
!!
!!if(totstr.gt.strn)goto 717
!!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!!write(9,*),'ifail',ifail,det_inv
!!print*,'Ifail',Ifail,i7
!!if(Ifail.eq.1)goto 717
!!!strno(totstr)=0
!!!totstr=totstr-1
!!!goto 317
!!!endif
!!!endif
!!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!!print*,'MatLDR_write_xmi',totstr
!!call MatLDR('str',strno,totstr,D)
!!
!!!print*,'DDD',(D(i),i=1,totstr)
!!ovlp=1.0
!!do i=1,totstr
!!ovlp=ovlp*D(i)
!!enddo
!!!print*,'ovval',1.0-ovlp,ovlp
!!if(1.0-ovlp.gt.ovval)goto 717
!!endif
!!
!!
!!m21=m21+nqset(m17)
!!if(m21.gt.strn) goto 717
!!!print*,'sourav_xmi_sym'
!!
!!!317 enddo
!!
!!!417 if(jj.eq.nqset(m17))then
!!
!!do m19=strset(m17-1)+1,strset(m17)
!!if(q_fac2(m19).ne.qul(m17))goto 617
!!i7=i7+1
!!print*,'new struc',i7,m19
!!do i3=1,nae
!!str2(i7,i3)=str3(m19,i3)
!!enddo
!!col(i7)=m19
!!qq(i7)=q_fac2(m19)
!!qq1(i7)=str_quality_1(m19)
!!qq2(i7)=str_quality_2(m19)
!!bondq4(i7)=bondq(m19)
!!flg=1
!!617 enddo
!!!print*,'i7 16',i7,m17
!!if(i7.eq.strn) then
!!!if(nfset.gt.0)then
!!ttqlty=0
!!tqlty=0
!!bqlty=0
!!sqlty=0
!!do m119=1,i7
!!ttqlty=ttqlty+qq(m119)
!!!tqlty=tqlty+qq1(m119)
!!!bqlty=bqlty+bondq4(m119)
!!!sqlty=sqlty+qq2(m119)
!!enddo
!!!endif
!!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!!if(ttqlty.le.ttqlty0)then
!!set_number=set_number+1
!!mns=mns+1
!!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!!if(nfset.eq.1.and.set_number.gt.1)then
!!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!!endif
!!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!enddo
!!!if(nfset.eq.0.or.nfset.gt.1) then
!!if(nfset.eq.3)then
!!Rumwrite=1
!!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!!endif
!!if(ovopt.eq.1) then
!!call MatLDR('str',col,i7,D)
!!
!!ovlp=1.0
!!do i=1,i7
!!ovlp=ovlp*D(i)
!!enddo
!!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!!endif
!!!write(9,910)'quality value',ttqlty
!!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!bqlty=0
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!!endif
!!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!!'ovval',1.0-ovlp
!!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!!if(mns.eq.max_set)then
!!if(u1.eq.mset-1) goto 379
!!close(9+u1)
!!mns=0
!!u1=u1+1
!!if(u1.le.9)then
!!write(a,'(I1)')u1
!!endif
!!if(u1.gt.9.and.u1.lt.100)then
!!write(a,'(I2)')u1
!!endif
!!if(u1.gt.99.and.u1.lt.1000)then
!!write(a,'(I3)')u1
!!endif
!!if(u1.eq.1000)then
!!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!!sets. Please increase the limit'
!!goto 379
!!endif
!!outfile='structure_set_'//trim(a)//trim('.dat')
!!open(unit=9+u1,file=outfile,status='unknown')
!!endif
!!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!!write(9,911),tqlty,sqlty,bqlty
!!!bqlty=0
!!!!set_number=set_number+1
!!!write(9,*),'Set_number=',set_number
!!!write(9,*),'    '
!!!endif
!!!endif
!!!if(nfset.eq.1)goto 379
!!endif
!!
!!717 enddo
!!!print*,'count',count
!!if(indpnt.eq.2)then
!!if (count.eq.1)goto 702
!!if (count.eq.2)goto 703
!!if (count.eq.3)goto 704
!!if (count.eq.4)goto 705
!!if (count.eq.5)goto 706
!!if (count.eq.6)goto 707
!!if (count.eq.7)goto 708
!!if (count.eq.8)goto 709
!!if (count.eq.9)goto 710
!!if (count.eq.10)goto 711
!!if (count.eq.11)goto 712
!!if (count.eq.12)goto 713
!!if (count.eq.13)goto 714
!!if (count.eq.14)goto 715
!!if (count.eq.15)goto 716
!!endif
!!!if (count.eq.16)goto 717
!!
!!if(indpnt.eq.1)then
!!if (count.eq.1)goto 701
!!if (count.eq.2)goto 702
!!if (count.eq.3)goto 703
!!if (count.eq.4)goto 704
!!if (count.eq.5)goto 705
!!if (count.eq.6)goto 706
!!if (count.eq.7)goto 707
!!if (count.eq.8)goto 708
!!if (count.eq.9)goto 709
!!if (count.eq.10)goto 710
!!if (count.eq.11)goto 711
!!if (count.eq.12)goto 712
!!if (count.eq.13)goto 713
!!if (count.eq.14)goto 714
!!if (count.eq.15)goto 715
!!if (count.eq.16)goto 716
!!endif
!!
!!716 enddo
!!
!!715 enddo
!!
!!714 enddo
!!
!!713 enddo
!!
!!712 enddo
!!
!!711 enddo
!!
!!710 enddo
!!
!!709 enddo
!!
!!708 enddo
!!
!!707 enddo
!!
!!706 enddo
!!
!!705 enddo
!!
!!704 enddo
!!
!!703 enddo
!!
!!702 enddo
!!
!!701 enddo
!!!enddo
!!
!!
!!!379 do m19=1,i7
!!!if(niao.eq.0)then
!!!write(9,900),qq(m19),(str2(m19,m20),m20=1,nae)
!!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!!write(9,901),qq(m19),bondq4(m19),1,':',niao,(str2(m19,m20),m20=1,nae)
!!!!write(9,901),m19,1,':',niao,(str2(m19,m20),m20=1,nae)
!!!!write(*,900),m19,niao,(str2(m19,m20),m20=1,nae)
!!!!write(9,*)qq(m19)
!!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!!write(9,909),qq(m19),1,':',niao,(str2(m19,m20),m20=1,nae)
!!!endif
!!
!!
!!!do m19=1,i7
!!!if(niao.eq.0)then
!!!write(9,900),qq1(m19),qq2(m19),bondq4(m19),'|',(str2(m19,m20),m20=1,nae)
!!!endif
!!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!!if(niao.gt.1)then
!!!write(9,901),qq1(m19),qq2(m19),bondq4(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!!endif
!!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!!if(niao.eq.1)then
!!!write(9,909),qq1(m19),qq2(m19),bondq4(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!!endif
!!!tqlty=tqlty+qq1(m19)
!!!bqlty=bqlty+bondq4(m19)
!!!sqlty=sqlty+qq2(m19)
!!!enddo
!!
!!!900 format(25I4)
!!!901 format(I3,x,I3,x,I1,a,I1,x,25I4)
!!!909 format(I3,x,I1,a,I1,x,25I4)
!!
!!900 format(I3,x,I3,x,I3,x,I3,x,a,x,25I4)
!!901 format(I3,x,I3,x,I3,x,I3,x,a,x,I1,a,I3,x,25I4)
!!909 format(I3,x,I3,x,I3,x,I3,x,a,x,I3,I3,x,25I4)
!!!910 format(a,x,a,x,a,x,a)
!!!911 format(15x,I3,7x,I3,7x,I3)
!!910 format(a,I3)
!!920 format(a,2x,I5,2x,a,2x,100I5)
!!911 format(15x,I3,7x,I3,7x,I3)
!!912 format(a,3x,F10.3)
!!
!!close(21)
!!
!!open(unit=121,file='script10',status='unknown')
!!write(121,111)'rm -rf','Rumer_Sets.dat'
!!close(121)
!!CALL SYSTEM ("chmod +x script10 ")
!!CALL SYSTEM ("./script10 ")
!!CALL SYSTEM ("rm script10 ")
!!111 format(a,x,a)
!!
!!231 format(50I3)
!!!deallocate(str2)
!!!deallocate(qq)
!!!deallocate(qq1)
!!!deallocate(qq2)
!!print*,'exit write_symm_xmi'
!!
!!379 return
!!end subroutine w_s_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_str_bond(nl,strn,str3,ncqs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! this subroutine creates the bond serials need to check the 
!! independency in bond slab way
use commondat
implicit none

integer::i,i1,i2,i3,j,j1,j2,nl,strn,ncqs,str3(15000,20),num_orb(20),num_orb1(20),&
totbnd,totseorb,str_bnd(15000,1000),str_rad(500,20),bond_sl(10000,2),least
common /chek/totbnd,totseorb,str_bnd,str_rad,bond_sl,num_orb


print*,'enter check_str_bond'
print*,'sourav1'
i1=0
do i=1,nae-nl*2
i1=i1+1
num_orb(i1)=str3(1,i+nl*2)
enddo

least=1000
do i2=1,i1
if(least.gt.num_orb(i2))least=num_orb(i2)
enddo

i3=0
320 do i2=1,i1
if(least.eq.num_orb(i2))then
i3=i3+1
num_orb1(i3)=least
goto 321
endif
enddo
321 least=least+1
if (i3.lt.i1)goto 320

do i=1,i1
num_orb(i)=num_orb1(i)
enddo

i2=0
do j1=1,i1
do j2=j1+1,i1
i2=i2+1
bond_sl(i2,1)=num_orb(j1)
bond_sl(i2,2)=num_orb(j2)
enddo
enddo

print*,'sourav2'
totbnd=i2
totseorb=i1

do i=1,15000
do i1=1,1000
str_bnd(i,i1)=0
enddo
enddo

do i=1,ncqs
do i1=nl*2+1,nae-nlast,2
do j=1,totbnd
j1=0
do i2=i1,i1+1
do i3=1,2
if(str3(i,i2).eq.bond_sl(j,i3))then
j1=j1+1
endif
enddo
enddo
if(j1.eq.2)then
str_bnd(i,j)=1
goto 200
endif
enddo
200 enddo
enddo

do i=1,500
do i1=1,20
str_rad(i,i2)=0
enddo
enddo

if(nlast.ne.0)then
do i=1,ncqs
do i1=nae-nlast+1,nae
do i2=1,totseorb
if(str3(i,i1).eq.num_orb(i2))then
str_rad(i,i2)=1
goto 201
endif
enddo
201 enddo
enddo
endif

!print*,'str',(str3(1,i),i=1,nae)
!do i=1,totbnd
!print*,'bond_sl',i,'<<',(bond_sl(i,j),j=1,2)
!enddo
do i=1,ncqs
write(*,100),'str3',i,'>',(str3(i,j),j=1,nae)
write(*,100),'str_bnd',i,'>',(str_bnd(i,j),j=1,totbnd)
enddo
100 format(a,I3,a,50I4)
!stop

print*,'exit check_str_bond'
return
end subroutine check_str_bond
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_xmi_new_2(nl,strn,str3,ncqs,q_fac2,quality_fac,ffvec2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig
common/str/str5,nstr7
common/qfacRumid/qq1,bondq4,qq2,qq,rumset

integer::nl,strn,ncqs,tostr,initstr,s,i,i1,i2,i3,i4,i5,i6,i7,i8,i9,m119,tqlty,set_number,&
mn1,mn2,&
m16,m15,m14,m13,m12,m11,m10,m9,m8,m7,m6,m5,m4,m3,m2,m1,m17,m18,m19,mm19,m20,m21,m23,m24,Ifailn&
,qul(100),nqul,j,jj,fg,flg,olp,qsn,countm19,i5up1,i5up2,i5up3,i5up4,i5up5,i5up6,i5up7,i5up8,i5up9 &
,i5up10,i5up11,i5up12,i5up13,i5up14,i5up15,i5up16,i5up17,i5up18,i5up19,i5up20,&
i1i7,i2i7,i3i7,i4i7,i5i7,i6i7,i7i7,i8i7,i9i7,i10i7,i11i7,i12i7,i13i7,i14i7,i15i7,i16i7,i17i7,i18i7,i19i7,i20i7&
,m1m21,m2m21,m3m21,m4m21,m5m21,m6m21,m7m21,m8m21,m9m21,m10m21,m11m21,m12m21,m13m21,m14m21,m15m21,m16m21&
,m17m21,m18m21,m19m21,m20m21,totstr,ifail,bqlty,sqlty,ttqlty1,ttqlty0,ttqlty2,ttqlty3,ttqlty,hqlty
integer::n19,l19,o19,j17,j27,j37,j47,m,mm121,mm221,mm321,mm421,x,y,det_inv,detc(100),detc1(100),set_num(100),u1,max_set,mns
integer::str3(15000,20),q_fac2(15000),finalvec(15000),strno(15000),sigsym(15000),tnqs_sig&
,ffvec2(15000,1000),bondq(15000),overlap(15000),qstrn(50),qqstrn(50),vec(15000,1000)&
,tndet,quality_fac(15000),str_quality_1(15000),str_quality_2(15000),tnqs,nssym,qulsym(15000)&
,symq(15000),col(1000),str2(2000,20),fvec(15000,1000),avg,bondc(100,50),numb,ct,bdpstr(100),sum1&
,bdset,str5(2000,20),nstr7,Rid,rumer(15000),rumer_rad(15000),iter_counter,nnn,nbd,bdet(50)&
,pref_radical(15000),mbondq(15000)
real*8::ovlp,std,var
real::bnd_std
Double Precision::D(1000),bond_det
character(len=10)::inact
character(len=6)::a
character(len=100)::outfile

integer::qq1(5000),qq2(5000),qq(5000),bondq4(5000),rumset
print*,'enter write_xmi_new_2',niao

mns=0
u1=1
max_set=75000

write(35,'(a)')'********** All the structures arranged serialy according to their qualities ********'
write(35,'(a)')'************************************************************************************'
if(nl.ne.0)write(35,122)'Structures in lone pair set ',nl,'is given below'
122 format(a,1x,I3,1x,a)
write(35,123)'sl no.','srtuctures','IAB','NN','SB','UDB','UDR','Tot'
if(radical.ne.1.and.prad.ne.0.and.nlast.ne.0.or.flg1.eq.1) call prio_rad_str(nl,str3,ncqs,pref_radical)
if(mnbond.ne.1.and.imbd.ne.0.or.flg1.eq.1) call main_bond_cal(nl,str3,ncqs,mbondq)
if(niao.eq.0)inact=''
if(niao.eq.1)inact='1  1'
if(niao.gt.1)write (inact, "(A2,I2)") "1:", niao
do i=1,ncqs
123 format(a,6x,a,8x,a,2x,a,2x,a,2x,a,2x,a,2x,a)
write(35,124),i,')',trim(inact),(str3(i,j),j=1,nae),str_quality_1(i),bondq(i),str_quality_2(i),mbondq(i),pref_radical(i),q_fac2(i)
enddo
124 format(I2,1x,a,2x,a,1x,30I3)
nbd=0
do i=1,50
bdet(i)=0
enddo

!!! 'check_str_bond' is necessary if we are checking independency with the bond
!!! slabe
 if(flg_ion.eq.0.and.vacorb.le.0)call check_str_bond(nl,strn,str3,ncqs)

iter_counter=0
if(nfset.eq.3.or.nfset.eq.5)then
rumset=0
call rumer_structures(nl,str3,ncqs,rumer,rumer_rad)
call write_rumer_xmi(nl,str3,ncqs,rumer,rumer_rad,quality_fac)
endif
ct=0
numb=0
std=0.0
do i=1,100
detc(i)=0
detc1(i)=0
enddo
avg=0
var=0.0

!do i=1,ncqs
!print*,'q_fac2:write_xmi_new_2',q_fac2(i)
!enddo

if(noq0.gt.strn)then
ttqlty0=noq0
else
ttqlty0=strn+noq0
endif
ttqlty1=strn+noq1
ttqlty2=strn+noq2
ttqlty3=strn+noq3

!print*,'nfset',nfset,ovopt
write(*,*)'ttqlty1,ttqlty2,ttqlty3',ttqlty1,ttqlty2,ttqlty3,ttqlty0,noq1,noq2,noq3,noq0,nfset
!
!print*,'ttqlty.le.ttqlty0',ttqlty,ttqlty0

!ovlpval=0.985

set_number=0
ttqlty=0
bqlty=0
tqlty=0
sqlty=0

jj=1
do m19=1,ncqs
if(m19.eq.1)qul(1)=q_fac2(1)
j=jj
do i=1,j
if(qul(i).eq.q_fac2(m19))goto 373
enddo
jj=jj+1
qul(i)=q_fac2(m19)
!print*,qul(i)
373 enddo
nqul=jj

do i=1,nqul
qsn=0
do m19=1,ncqs
if(qul(i).ne.q_fac2(m19))goto 374
qsn=qsn+1
374 enddo
qqstrn(i)=qsn
enddo
qstrn(1)=0
j=1
do i=1,nqul
j=j+1
s=qqstrn(i)
qstrn(j)=qstrn(j-1)+s
enddo


!print*,'nqul',nqul
!print*,'last_sourav'
do m19=1,ncqs
!!!!print*,m19
write(*,231)(str3(m19,i4),i4=1,nae)
!print*,'q_fac2',q_fac2(m19)
!print*,'bondq',bondq(m19)
!print*,'qul',qul(m19)
!!print*,'niao',niao,strn
!!write(*,231)(ffvec2(m19,i4),i4=1,ndet)
!!write(*,231)(bond_count(m19,i4),i4=1,nmbond)
!!!!write(*,231)(ffvec2(m19,i4),i4=1,ndet)
enddo
!uosn=uoptstr
!endif


434 format (40I3)

do m19=1,ncqs

!if(max_set.eq.mns)then
!print*,'mns',mns
!u1=u1+1
!write(a,'(I4)')u1
!outfile=trim('structure_set_')//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif

!!!! quality 1 structures !!!!!!!
!write(*,*)'111111111111111111111',m19,ncqs
totstr=0
i7=0
m21=0

do m20=1,10000
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo



totstr=totstr+1
strno(totstr)=m19
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifail)
if(Ifail.eq.1)then
!iter_counter=iter_counter+1
!print*,'iter_counter',iter_counter,Ifail
goto 200
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(totstr.gt.1.and.ovopt.eq.vpt) then
call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 200
endif
m21=m21+1
if(m21.gt.strn) goto 316
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
!print*,'i7',i7,m19
!print*,'sourav1',q_fac2(m19)
qq(i7)=q_fac2(m19)
!print*,'sourav2'
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)

!print*,'nmbond',nmbond
!do i3=1,nmbond
!bondc(i7,i3)=bond_count(m19,i3)
!enddo
!print*,'bondq4(i7)',bondq4(i7),m19

if(i7.eq.strn) then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
endif

do m119=1,i7
if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(bondc(m119,i4),i4=1,nmbond)
!do i4=1,nmbond
!detc(i4)=detc(i4)+bondc(m119,i4)
!enddo
!do i4=1,nmbond
!bdpstr(m119)=bdpstr(m119)+bondc(m119,i4)
!enddo

enddo
if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif
if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif
if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
goto 200
endif


i1i7=i7
m1m21=m21


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!loop 2 start
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do m18=m19+1,ncqs
!write(*,*)'222222222222222222222222222',m18,ncqs

do m20=i1i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo
i7=i1i7
totstr=i1i7
m21=m1m21

totstr=totstr+1
strno(totstr)=m18
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifail)
if(Ifail.eq.1)then
goto 201
endif
!!!!!!!!!!!!!!!!!!
if(totstr.gt.1.and.ovopt.eq.vpt) then
call MatLDR('str',strno,totstr,D)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 201
endif
m21=m21+1
if(m21.gt.strn) goto 316

i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m18,i3)
enddo
col(i7)=m18
!print*,'i7',i7,m18
qq(i7)=q_fac2(m18)
qq1(i7)=str_quality_1(m18)
qq2(i7)=str_quality_2(m18)
bondq4(i7)=bondq(m18)

!do i3=1,nmbond
!bondc(i7,i3)=bond_count(m18,i3)
!enddo

if(i7.eq.strn) then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
endif

do i4=1,100
detc(i4)=0
detc1(i4)=0
bdpstr(i4)=0
enddo
do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(bondc(m119,i4),i4=1,nmbond)
!do i4=1,nmbond
!detc(i4)=detc(i4)+bondc(m119,i4)
!enddo
!do i4=1,nmbond
!bdpstr(m119)=bdpstr(m119)+bondc(m119,i4)
!enddo


enddo


if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif
if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!bdset=0
!do i=1,i7
!if(bdpstr(i).ne.0)then
!bdset=bdset+1
!endif
!enddo
!numb=nmbond
!write(9,232)'bond',(detc(i),i=1,numb)
!ct=0
!do i=1,numb
!ct=ct+1
!detc1(ct)=detc(i)
!enddo
!numb=ct
!sum1=0
!do i=1,numb
!avg=avg+detc1(i)
!sum1=sum1+detc1(i)
!enddo
!avg=avg/(numb)
!do i=1,numb
!var=var+(detc1(i)-avg)**2
!enddo
!var=var/(numb)
!std=sqrt(var)
!write(9,233)'STD_BOND',',',set_number,',',std,',',(detc(i),i=1,numb),',',sum1,',',(bdpstr(i),i=1,i7)&
!,',',bdset
!233 format(a,a,I5,a,F10.3,a,6I3,a,I3,a,2I3,a,I3)
!avg=0
!var=0.0

bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif


!if(nfset.eq.1) write(9,*),'ovval',1.0-ovlp
!if(nfset.eq.1) goto 316
!if(m18.eq.ncqs)goto 200

goto 201
endif



i5up2=i5
i2i7=i7
m2m21=m21


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! loop 3 start
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do m17=m18+1,ncqs

!!!! quality 3 structures !!!!!!!
!write(*,*)'333333333333333333333333333',m17,ncqs

do m20=i2i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo
i7=i2i7
totstr=i2i7
m21=m2m21

totstr=totstr+1
strno(totstr)=m17
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
!call check_ind(str3,totstr,nl,ncqs,strno,Ifail)

if(Ifail.eq.1)then

!iter_counter=iter_counter+1
!print*,'iter_counter',iter_counter,Ifail
!call bond_det_vec(strno,totstr,bond_det)

goto 202
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(totstr.gt.1.and.ovopt.eq.vpt) then
call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo

!print*,'ovval',1.0-ovlp,ovlp

if(1.0-ovlp.gt.ovval)goto 202
endif
m21=m21+1
if(m21.gt.strn) goto 316

i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m17,i3)
enddo
col(i7)=m17
!print*,'i7',i7,m17
qq(i7)=q_fac2(m17)
qq1(i7)=str_quality_1(m17)
qq2(i7)=str_quality_2(m17)
bondq4(i7)=bondq(m17)

!do i3=1,nmbond
!fvec(i7,i3)=ffvec2(m17,i3)
!bondc(i7,i3)=bond_count(m17,i3)
!enddo

if(i7.eq.strn) then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
endif

do m119=1,i7


if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

enddo

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif

!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then
!print*,'ovopt',ovopt

if(ovopt.eq.1) then

!print*,'MatLDR3'
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif
if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

!if(m17.eq.ncqs)goto 201

goto 202
endif


i5up3=i5
i3i7=i7
m3m21=m21

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! loop 4 start
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do m16=m17+1,ncqs

!!!! quality 4 structures !!!!!!!
!write(*,*)'4444444444444444444444444',m16,ncqs


do m20=i3i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i3i7
totstr=i3i7
m21=m3m21



totstr=totstr+1
strno(totstr)=m16
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn

if(Ifail.eq.1)then

!iter_counter=iter_counter+1
!print*,'iter_counter',iter_counter,Ifail
!call bond_det_vec(strno,totstr,bond_det)

goto 203
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(totstr.gt.1.and.ovopt.eq.vpt) then
call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 203
endif
m21=m21+1
if(m21.gt.strn) goto 316

i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m16,i3)
enddo
col(i7)=m16

!print*,'i7',i7,m16

qq(i7)=q_fac2(m16)
qq1(i7)=str_quality_1(m16)
qq2(i7)=str_quality_2(m16)
bondq4(i7)=bondq(m16)

!do i3=1,nmbond
!!fvec(i7,i3)=ffvec2(m16,i3)
!bondc(i7,i3)=bond_count(m16,i3)
!enddo

if(i7.eq.strn)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
endif

do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(fvec(m119,i4),i4=1,ndet)
!do i4=1,ndet
!detc(fvec(m119,i4))=detc(fvec(m119,i4))+1
!enddo



enddo

if(flgst.eq.1.or.flgst.eq.3) then
call ion_struc
endif
!if(nfset.eq.3)call Rumer_set_id(str2,i7,nl,Rid,set_num)
if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif

!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then
!print*,'ovopt',ovopt

if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!write(9,232)'det',(detc(i),i=1,tot_ndet-1)
!do i=1,tot_ndet-1
!avg=avg+detc(i)
!enddo
!avg=avg/(tot_ndet-1)
!do i=1,tot_ndet-1
!var=var+(detc(i)-avg)**2
!enddo
!var=var/(tot_ndet-1)
!std=sqrt(var)
!write(9,233)'STD_DET',std
!do i=1,1000
!detc(i)=0
!enddo
!avg=0
!var=0.0
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0

write(9,*),'Set_number=',set_number
write(9,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

!if(nfset.eq.1) write(9,*),'ovval',1.0-ovlp
!if(nfset.eq.1) goto 316
!if(m16.eq.ncqs)goto 202

goto 203
endif


!152 enddo

!if(noq.eq.4)goto 203

i4i7=i7
m4m21=m21

do m15=m16+1,ncqs

!!! quality 5 structures !!!!!!!
!write(*,*)'55555555555555555555555555',m15,ncqs
!iter_counter=iter_counter+1
!print*,'iter_counter',iter_counter
!if(max_set.eq.mns)then
!print*,'mns',mns
!u1=u1+1
!write(a,'(I4)')u1
!outfile=trim('structure_set_')//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif

do m20=i4i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i4i7
totstr=i4i7
m21=m4m21



!do m19=m14,ncqs


!if(q_fac2(m19).ne.q_fac2(m14))goto 142

totstr=totstr+1
strno(totstr)=m15
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn

if(Ifail.eq.1)then
!!iter_counter=iter_counter+1
goto 204
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if(nfset.eq.1.and.totstr.gt.1) then

if(totstr.gt.1.and.ovopt.eq.vpt) then

!print*,'MatLDR_write_xmi',totstr
!print*,'MatLDR1'

call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 204
endif
m21=m21+1
if(m21.gt.strn) goto 316

i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m15,i3)
enddo
col(i7)=m15
!print*,'i7',i7,m15

qq(i7)=q_fac2(m15)
qq1(i7)=str_quality_1(m15)
qq2(i7)=str_quality_2(m15)
bondq4(i7)=bondq(m15)

!do i3=1,nmbond
!bondc(i7,i3)=bond_count(m15,i3)
!!fvec(i7,i3)=ffvec2(m15,i3)
!enddo

if(i7.eq.strn)then

!write(9,*),'srt',(strno(m119),m119=1,i7)
!ttqlty=0
!if(nfset.gt.0)then

ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo

!endif
!if(set_number.eq.1)hqlty=ttqlty
!if(set_number.eq.1.and.nfset.eq.3)hqlty=ttqlty+noq
!write(*,*)'hqlty,ttqlty',hqlty,ttqlty
!if(ttqlty.le.hqlty)then
!if(ttqlty1.ge.tqlty.and.ttqlty2.ge.bqlty.and.ttqlty3.ge.sqlty.or.nfset.eq.0)then
!print*,'tqlty,bqlty,sqlty,ttqlty',tqlty,bqlty,sqlty,ttqlty
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
mns=mns+1
iter_counter=iter_counter+1
print*,'iter_counter',iter_counter,Ifail
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do i4=1,100
detc(i4)=0
detc1(i4)=0
bdpstr(i4)=0
enddo


do m119=1,i7

if(niao.eq.0)then
write(9+u1,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
!write(*,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(bondc(m119,i4),i4=1,nmbond)
!do i4=1,nmbond
!detc(i4)=detc(i4)+bondc(m119,i4)
!enddo
!do i4=1,nmbond
!bdpstr(m119)=bdpstr(m119)+bondc(m119,i4)
!enddo
!print*,'bond',(detc(i4),i4=1,nmbond)

enddo
!stop


call bond_det_vec(nl,i7,str2,bnd_std)

!!do i=1,nbd
!!if(bdet(i).eq.int(bond_det))goto 777
!!enddo
!!nbd=nbd+1
!!bdet(nbd)=int(bond_det)
!!777 write(9,913)'bond_vec_det =',bond_det,Ifail,Ifailn
!write(9,*),'Ifail_test_check_ind',Ifail,Ifailn
!if(nfset.eq.3)call Rumer_set_id(str2,i7,nl,Rid,set_num)

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif

!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then

if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9+u1,921)'quality value',ttqlty,'bnd_std',bnd_std
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!numb=tot_ndet-1
!bdset=0
!do i=1,i7
!if(bdpstr(i).ne.0)then
!bdset=bdset+1
!endif
!enddo
!!numb=nmbond
!!write(9,232)'bond',(detc(i),i=1,numb)
!ct=0
!do i=1,numb
!!if(detc(i).eq.0)goto 333
!ct=ct+1
!detc1(ct)=detc(i)
!!endif
!333 enddo
!numb=ct
!!print*,'ct',ct,(detc1(i),i=1,ct)
!sum1=0
!do i=1,numb
!avg=avg+detc1(i)
!sum1=sum1+detc1(i)
!enddo
!avg=avg/(numb)
!do i=1,numb
!var=var+(detc1(i)-avg)**2
!enddo
!var=var/(numb)
!std=sqrt(var)
!!write(9,233)'STD_BOND',',',set_number,',',std,',',(detc(i),i=1,numb),',',sum1,',',(bdpstr(i),i=1,i7)&
!!,',',bdset
!avg=0
!var=0.0
!stop
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0

!set_number=set_number+1
!write(9,*),'Set_number=',set_number,iter_counter,Ifail

write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

!if(nfset.eq.1) write(9,*),'ovlpval',1.0-ovlp
!if(nfset.eq.1) goto 316
!if(m15.eq.ncqs)goto 203

if(mns.eq.max_set)then
!print*,'mns',mns
!print*,'u1',u1
close(9+u1)
mns=0
u1=u1+1
!write(a,'(I1)')u1
!a=trim(a)
!print*,'a',a
outfile='structure_set_'//trim(int_num(u1))//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

goto 204
endif

!142 enddo
!if(noq.eq.5)goto 204

i5i7=i7
m5m21=m21

!do m13=1+qstrn(6),qstrn(7)
do m14=m15+1,ncqs

!!!! quality 6 structures !!!!!!!
!write(*,*)'66666666666666666666666666'


do m20=i5i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i5i7
totstr=i5i7
m21=m5m21



!do m19=m13,ncqs
!if(q_fac2(m19).ne.q_fac2(m13))goto 132

totstr=totstr+1
strno(totstr)=m14
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail6',ifail

if(Ifail.eq.1)then
!strno(totstr)=0
!totstr=totstr-1
goto 205
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if(nfset.eq.1.and.totstr.gt.1) then

if(totstr.gt.1.and.ovopt.eq.vpt) then

!print*,'MatLDR_write_xmi',totstr
!print*,'MatLDR1'

call MatLDR('str',strno,totstr,D)
!print*,'DDD',(D(i),i=1,totstr)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 205
endif
m21=m21+1

if(m21.gt.strn) goto 316

i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m14,i3)
enddo
col(i7)=m14

!print*,'i7',i7,m14

qq(i7)=q_fac2(m14)
qq1(i7)=str_quality_1(m14)
qq2(i7)=str_quality_2(m14)
bondq4(i7)=bondq(m14)

!do i3=1,nmbond
!bondc(i7,i3)=bond_count(m14,i3)
!enddo
!do i3=1,ndet
!fvec(i7,i3)=ffvec2(m14,i3)
!enddo

if(i7.eq.strn)then

!write(9,*),'srt',(strno(m119),m119=1,i7)
!ttqlty=0
!if(nfset.gt.0)then

ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo

!endif
!if(set_number.eq.1)hqlty=ttqlty
!if(set_number.eq.1.and.nfset.eq.3)hqlty=ttqlty+noq
!write(*,*)'hqlty,ttqlty',hqlty,ttqlty
!if(ttqlty.le.hqlty)then
!if(ttqlty1.ge.tqlty.and.ttqlty2.ge.bqlty.and.ttqlty3.ge.sqlty.or.nfset.eq.0)then
!print*,'tqlty,bqlty,sqlty,ttqlty',tqlty,bqlty,sqlty,ttqlty
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty

!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty

!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

endif

do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(fvec(m119,i4),i4=1,ndet)
!do i4=1,ndet
!detc(fvec(m119,i4))=detc(fvec(m119,i4))+1
!enddo


enddo

!if(nfset.eq.3)call Rumer_set_id(str2,i7,nl,Rid,set_num)

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif

!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then

if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!write(9,232)'det',(detc(i),i=1,tot_ndet-1)
!do i=1,tot_ndet-1
!avg=avg+detc(i)
!enddo
!avg=avg/(tot_ndet-1)
!do i=1,tot_ndet-1
!var=var+(detc(i)-avg)**2
!enddo
!var=var/(tot_ndet-1)
!std=sqrt(var)
!write(9,233)'STD_DET',std
!do i=1,1000
!detc(i)=0
!enddo
!avg=0
!var=0.0
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

!if(nfset.eq.1) write(9,*),'ovlpval',1.0-ovlp
!if(nfset.eq.1) goto 316

goto 205
endif

!132 enddo

!if(noq.eq.6)goto 205

i6i7=i7
m6m21=m21

do m13=m14+1,ncqs

!!!! quality 7 structures !!!!!!!
!write(*,*)'7777777777777777777777777777'


do m20=i6i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i6i7
totstr=i6i7
m21=m6m21


!if(q_fac2(m19).ne.q_fac2(m12))goto 122

totstr=totstr+1
strno(totstr)=m13
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail7',ifail

if(Ifail.eq.1)then
!strno(totstr)=0
!totstr=totstr-1
goto 206
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if(nfset.eq.1.and.totstr.gt.1) then

if(totstr.gt.1.and.ovopt.eq.vpt) then

!print*,'MatLDR_write_xmi',totstr
!print*,'MatLDR1'

call MatLDR('str',strno,totstr,D)
!print*,'DDD',(D(i),i=1,totstr)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 206
endif
m21=m21+1

if(m21.gt.strn) goto 316


i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m13,i3)
enddo
col(i7)=m13
!print*,'i7',i7,m13

qq(i7)=q_fac2(m13)
qq1(i7)=str_quality_1(m13)
qq2(i7)=str_quality_2(m13)
bondq4(i7)=bondq(m13)

!do i3=1,nmbond
!bondc(i7,i3)=bond_count(m13,i3)
!enddo
!do i3=1,ndet
!fvec(i7,i3)=ffvec2(m13,i3)
!enddo

if(i7.eq.strn)then

!write(9,*),'srt',(strno(m119),m119=1,i7)
!ttqlty=0
!if(nfset.gt.0)then

ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo

!endif
!if(set_number.eq.1)hqlty=ttqlty
!if(set_number.eq.1.and.nfset.eq.3)hqlty=ttqlty+noq
!write(*,*)'hqlty,ttqlty',hqlty,ttqlty
!if(ttqlty.le.hqlty)then
!if(ttqlty1.ge.tqlty.and.ttqlty2.ge.bqlty.and.ttqlty3.ge.sqlty.or.nfset.eq.0)then
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(fvec(m119,i4),i4=1,ndet)
!do i4=1,ndet
!detc(fvec(m119,i4))=detc(fvec(m119,i4))+1
!enddo


enddo

!if(nfset.eq.3)call Rumer_set_id(str2,i7,nl,Rid,set_num)

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif

!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then

if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!write(9,232)'det',(detc(i),i=1,tot_ndet-1)
!do i=1,tot_ndet-1
!avg=avg+detc(i)
!enddo
!avg=avg/(tot_ndet-1)
!do i=1,tot_ndet-1
!var=var+(detc(i)-avg)**2
!enddo
!var=var/(tot_ndet-1)
!std=sqrt(var)
!write(9,233)'STD_DET',std
!do i=1,1000
!detc(i)=0
!enddo
!avg=0
!var=0.0
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0

write(9,*),'Set_number=',set_number
write(9,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

!if(nfset.eq.1) write(9,*),'ovlpval',1.0-ovlp
!if(nfset.eq.1) goto 316
goto 206
endif


!if(noq.eq.7)goto 206

i7i7=i7
m7m21=m21

do m12=m13+1,ncqs

!!!! quality 8 structures !!!!!!!
!write(*,*)'8888888888888888888888888888888'


do m20=i7i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i7i7
totstr=i7i7
m21=m7m21

!if(q_fac2(m19).ne.q_fac2(m11))goto 112

totstr=totstr+1
strno(totstr)=m12
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail8',ifail,(strno(i),i=1,totstr)

if(Ifail.eq.1)then
!strno(totstr)=0
!totstr=totstr-1
goto 207
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if(nfset.eq.1.and.totstr.gt.1) then

if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
!print*,'MatLDR1'
call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 207
endif
m21=m21+1

if(m21.gt.strn) goto 316

i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m12,i3)
enddo
col(i7)=m12
!print*,'i7',i7,m12

qq(i7)=q_fac2(m12)
qq1(i7)=str_quality_1(m12)
qq2(i7)=str_quality_2(m12)
bondq4(i7)=bondq(m12)

!do i3=1,nmbond
!bondc(i7,i3)=bond_count(m12,i3)
!enddo
!do i3=1,ndet
!fvec(i7,i3)=ffvec2(m12,i3)
!enddo

if(i7.eq.strn)then

!write(9,*),'srt',(strno(m119),m119=1,i7)
!ttqlty=0
!hqlty=0
!if(nfset.gt.0)then

ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo

!endif
!if(set_number.eq.1)hqlty=ttqlty
!if(set_number.eq.1.and.nfset.eq.3)hqlty=ttqlty+noq
!if(ttqlty.le.hqlty)then
!if(ttqlty1.ge.tqlty.and.ttqlty2.ge.bqlty.and.ttqlty3.ge.sqlty.or.nfset.eq.0)then
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(fvec(m119,i4),i4=1,ndet)
!do i4=1,ndet
!detc(fvec(m119,i4))=detc(fvec(m119,i4))+1
!enddo



enddo

!if(nfset.eq.3)call Rumer_set_id(str2,i7,nl,Rid,set_num)

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif

!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then

if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!write(9,232)'det',(detc(i),i=1,tot_ndet-1)
!do i=1,tot_ndet-1
!avg=avg+detc(i)
!enddo
!avg=avg/(tot_ndet-1)
!do i=1,tot_ndet-1
!var=var+(detc(i)-avg)**2
!enddo
!var=var/(tot_ndet-1)
!std=sqrt(var)
!write(9,233)'STD_DET',std
!do i=1,1000
!detc(i)=0
!enddo
!avg=0
!var=0.0
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

!if(nfset.eq.1) write(9,*),'ovlpval',1.0-ovlp
!if(nfset.eq.1) goto 316

goto 207
endif


!if(noq.eq.8)goto 207

i8i7=i7
m8m21=m21

do m11=m12+1,ncqs

!!!! quality 9 structures !!!!!!!
!write(*,*)'99999999999999999999999999'


do m20=i8i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i8i7
totstr=i8i7
m21=m8m21



!do m19=m10,ncqs
!if(q_fac2(m19).ne.q_fac2(m10))goto 102

totstr=totstr+1
strno(totstr)=m11
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail9',ifail,(strno(i),i=1,totstr)

if(Ifail.eq.1)then
!strno(totstr)=0
!totstr=totstr-1
goto 208
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
!print*,'MatLDR1'

call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 208
endif
m21=m21+1

if(m21.gt.strn) goto 316


i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m11,i3)
enddo
col(i7)=m11
!print*,'i7',i7,m11

qq(i7)=q_fac2(m11)
qq1(i7)=str_quality_1(m11)
qq2(i7)=str_quality_2(m11)
bondq4(i7)=bondq(m11)

!do i3=1,nmbond
!bondc(i7,i3)=bond_count(m11,i3)
!enddo
!do i3=1,ndet
!fvec(i7,i3)=ffvec2(m11,i3)
!enddo

if(i7.eq.strn)then

!write(9,*),'srt',(strno(m119),m119=1,i7)
!ttqlty=0
!hqlty=0
!if(nfset.gt.0)then

ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo

!endif
!if(set_number.eq.1)hqlty=ttqlty
!if(set_number.eq.1.and.nfset.eq.3)hqlty=ttqlty+noq
!if(ttqlty.le.hqlty)then
!if(ttqlty1.ge.tqlty.and.ttqlty2.ge.bqlty.and.ttqlty3.ge.sqlty.or.nfset.eq.0)then
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(fvec(m119,i4),i4=1,ndet)
!do i4=1,ndet
!detc(fvec(m119,i4))=detc(fvec(m119,i4))+1
!enddo



enddo

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif

!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then

if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!write(9,232)'det',(detc(i),i=1,tot_ndet-1)
!do i=1,tot_ndet-1
!avg=avg+detc(i)
!enddo
!avg=avg/(tot_ndet-1)
!do i=1,tot_ndet-1
!var=var+(detc(i)-avg)**2
!enddo
!var=var/(tot_ndet-1)
!std=sqrt(var)
!write(9,233)'STD_DET',std
!do i=1,1000
!detc(i)=0
!enddo
!avg=0
!var=0.0
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!if(nfset.eq.1) write(9,*),'ovval',1.0-ovlp
!if(nfset.eq.1) goto 316
goto 208
endif


i9i7=i7
m9m21=m21

do m10=m11+1,ncqs

!!!! quality 10 structures !!!!!!!
!write(*,*)'1010101010101010101010'


do m20=i9i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i9i7
totstr=i9i7
m21=m9m21

totstr=totstr+1
strno(totstr)=m10
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail10',ifail,(strno(i),i=1,totstr)

if(Ifail.eq.1)then
!strno(totstr)=0
!totstr=totstr-1
goto 209
endif

!if(nfset.eq.1.and.totstr.gt.1) then
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
!print*,'MatLDR1'
call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 209
endif
m21=m21+1

if(m21.gt.strn) goto 316


i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m10,i3)
enddo
col(i7)=m10
!print*,'i7',i7,m10

qq(i7)=q_fac2(m10)
qq1(i7)=str_quality_1(m10)
qq2(i7)=str_quality_2(m10)
bondq4(i7)=bondq(m10)

!do i3=1,nmbond
!bondc(i7,i3)=bond_count(m10,i3)
!enddo
!do i3=1,ndet
!fvec(i7,i3)=ffvec2(m10,i3)
!enddo

if(i7.eq.strn)then

!write(9,*),'srt',(strno(m119),m119=1,i7)
!ttqlty=0
!hqlty=0
!if(nfset.gt.0)then

ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo

!endif
!if(set_number.eq.1)hqlty=ttqlty
!if(set_number.eq.1.and.nfset.eq.3)hqlty=ttqlty+noq
!if(ttqlty.le.hqlty)then
!if(ttqlty1.ge.tqlty.and.ttqlty2.ge.bqlty.and.ttqlty3.ge.sqlty.or.nfset.eq.0)then
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(fvec(m119,i4),i4=1,ndet)
!do i4=1,ndet
!detc(fvec(m119,i4))=detc(fvec(m119,i4))+1
!enddo



enddo

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif

!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then

if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!write(9,232)'det',(detc(i),i=1,tot_ndet-1)
!do i=1,tot_ndet-1
!avg=avg+detc(i)
!enddo
!avg=avg/(tot_ndet-1)
!do i=1,tot_ndet-1
!var=var+(detc(i)-avg)**2
!enddo
!var=var/(tot_ndet-1)
!std=sqrt(var)
!write(9,233)'STD_DET',std
!do i=1,1000
!detc(i)=0
!enddo
!avg=0
!var=0.0
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

!if(nfset.eq.1) write(9,*),'ovval',1.0-ovlp
!if(nfset.eq.1) goto 316
goto 209
endif


i10i7=i7
m10m21=m21

do m9=m10+1,ncqs

!!!! quality 11 structures !!!!!!!
!write(*,*)'11 11 11 11 11 11 11'


do m20=i10i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i10i7
totstr=i10i7
m21=m10m21

totstr=totstr+1
strno(totstr)=m9
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail11',ifail,(strno(i),i=1,totstr)

if(Ifail.eq.1)then
!strno(totstr)=0
!totstr=totstr-1
goto 210
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
!print*,'MatLDR1'
call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 210
endif
m21=m21+1

if(m21.gt.strn) goto 316

i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m9,i3)
enddo
col(i7)=m9

!print*,'i7*****',i7,m9

qq(i7)=q_fac2(m9)
qq1(i7)=quality_fac(m9)
qq1(i7)=str_quality_1(m9)
qq2(i7)=str_quality_2(m9)
bondq4(i7)=bondq(m9)

!do i3=1,nmbond
!bondc(i7,i3)=bond_count(m9,i3)
!enddo
!do i3=1,ndet
!fvec(i7,i3)=ffvec2(m9,i3)
!enddo

if(i7.eq.strn)then

!write(9,*),'srt',(strno(m119),m119=1,i7)
!ttqlty=0
!hqlty=0
!if(nfset.gt.0)then

ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo

!if(set_number.eq.1)hqlty=ttqlty
!if(set_number.eq.1.and.nfset.eq.3)hqlty=ttqlty+noq
!if(ttqlty.le.hqlty)then
!if(ttqlty1.ge.tqlty.and.ttqlty2.ge.bqlty.and.ttqlty3.ge.sqlty.or.nfset.eq.0)then
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(fvec(m119,i4),i4=1,ndet)
!do i4=1,ndet
!detc(fvec(m119,i4))=detc(fvec(m119,i4))+1
!enddo

enddo

!strno(totstr)=0
!totstr=totstr-1
if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif

!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!write(9,232)'det',(detc(i),i=1,tot_ndet-1)
!do i=1,tot_ndet-1
!avg=avg+detc(i)
!enddo
!avg=avg/(tot_ndet-1)
!do i=1,tot_ndet-1
!var=var+(detc(i)-avg)**2
!enddo
!var=var/(tot_ndet-1)
!std=sqrt(var)
!write(9,233)'STD_DET',std
!do i=1,1000
!detc(i)=0
!enddo
!avg=0
!var=0.0
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

!if(nfset.eq.1) write(9,*),'ovval',1.0-ovlp
!if(nfset.eq.1) goto 316

goto 210
endif

i11i7=i7
m11m21=m21

do m8=m9+1,ncqs

!!!! quality 10 structures !!!!!!!
!write(*,*)'12 12 12 12 12 12'


do m20=i11i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i11i7
totstr=i11i7
m21=m11m21


totstr=totstr+1
strno(totstr)=m8
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail10',ifail,(strno(i),i=1,totstr)

if(Ifail.eq.1)then

!strno(totstr)=0
!totstr=totstr-1
goto 211
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
if(totstr.gt.1.and.ovopt.eq.vpt) then

!print*,'MatLDR_write_xmi',totstr
!print*,'MatLDR1'
call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 211
endif
m21=m21+1

if(m21.gt.strn) goto 316


i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m8,i3)
enddo
col(i7)=m8
!print*,'i7',i7,m8

qq(i7)=q_fac2(m8)
qq1(i7)=str_quality_1(m8)
qq2(i7)=str_quality_2(m8)
bondq4(i7)=bondq(m8)

!do i3=1,nmbond
!bondc(i7,i3)=bond_count(m8,i3)
!enddo
!do i3=1,ndet
!fvec(i7,i3)=ffvec2(m8,i3)
!enddo

if(i7.eq.strn)then

!write(9,*),'srt',(strno(m119),m119=1,i7)
!ttqlty=0
!hqlty=0
!if(nfset.gt.0)then

ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo

!endif
!if(set_number.eq.1)hqlty=ttqlty
!if(set_number.eq.1.and.nfset.eq.3)hqlty=ttqlty+noq
!!if(ttqlty.le.hqlty)then
!if(ttqlty1.ge.tqlty.and.ttqlty2.ge.bqlty.and.ttqlty3.ge.sqlty.or.nfset.eq.0)then
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
endif

do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(fvec(m119,i4),i4=1,ndet)
!do i4=1,ndet
!detc(fvec(m119,i4))=detc(fvec(m119,i4))+1
!enddo



enddo

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif

!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!write(9,232)'det',(detc(i),i=1,tot_ndet-1)
!do i=1,tot_ndet-1
!avg=avg+detc(i)
!enddo
!avg=avg/(tot_ndet-1)
!do i=1,tot_ndet-1
!var=var+(detc(i)-avg)**2
!enddo
!var=var/(tot_ndet-1)
!std=sqrt(var)
!write(9,233)'STD_DET',std
!do i=1,1000
!detc(i)=0
!enddo
!avg=0
!var=0.0
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

!if(nfset.eq.1) write(9,*),'ovval',1.0-ovlp
!if(nfset.eq.1) goto 316

goto 211
endif

i12i7=i7
m12m21=m21

do m7=m8+1,ncqs

!!!! quality 10 structures !!!!!!!
!write(*,*)'13 13 13 13 13 13'


do m20=i12i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i12i7
totstr=i12i7
m21=m12m21



totstr=totstr+1
strno(totstr)=m7
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail10',ifail,(strno(i),i=1,totstr)

if(Ifail.eq.1)then
!strno(totstr)=0
!totstr=totstr-1
goto 212
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
!print*,'MatLDR1'
call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 212
endif
m21=m21+1

if(m21.gt.strn) goto 316


i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m7,i3)
enddo
col(i7)=m7
!print*,'i7',i7,m7

qq(i7)=q_fac2(m7)
qq1(i7)=str_quality_1(m7)
qq2(i7)=str_quality_2(m7)
bondq4(i7)=bondq(m7)

!do i3=1,nmbond
!bondc(i7,i3)=bond_count(m7,i3)
!enddo
!do i3=1,ndet
!fvec(i7,i3)=ffvec2(m7,i3)
!enddo

if(i7.eq.strn)then

!write(9,*),'srt',(strno(m119),m119=1,i7)
!ttqlty=0
!hqlty=0

ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo

!if(set_number.eq.1)hqlty=ttqlty
!if(set_number.eq.1.and.nfset.eq.3)hqlty=ttqlty+noq
!if(ttqlty.le.hqlty)then
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!if(ttqlty1.ge.tqlty.and.ttqlty2.ge.bqlty.and.ttqlty3.ge.sqlty.or.nfset.eq.0)then

!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
endif

do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(fvec(m119,i4),i4=1,ndet)
!do i4=1,ndet
!detc(fvec(m119,i4))=detc(fvec(m119,i4))+1
!enddo



enddo

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif
!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!write(9,232)'det',(detc(i),i=1,tot_ndet-1)
!do i=1,tot_ndet-1
!avg=avg+detc(i)
!enddo
!avg=avg/(tot_ndet-1)
!do i=1,tot_ndet-1
!var=var+(detc(i)-avg)**2
!enddo
!var=var/(tot_ndet-1)
!std=sqrt(var)
!write(9,233)'STD_DET',std
!do i=1,1000
!detc(i)=0
!enddo
!avg=0
!var=0.0
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!if(nfset.eq.1) write(9,*),'ovval',1.0-ovlp
!if(nfset.eq.1) goto 316
goto 212
endif

i13i7=i7
m13m21=m21

do m6=m7+1,ncqs

!!!! quality 10 structures !!!!!!!
!write(*,*)'14 14 14 14 14 14'


do m20=i13i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i13i7
totstr=i13i7
m21=m13m21

totstr=totstr+1
strno(totstr)=m6
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail10',ifail,(strno(i),i=1,totstr)

if(Ifail.eq.1)then
!strno(totstr)=0
!totstr=totstr-1
goto 213
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
!print*,'MatLDR1'
call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 213
endif
m21=m21+1

if(m21.gt.strn) goto 316


i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m6,i3)
enddo
col(i7)=m6
!print*,'i7',i7,m6

qq(i7)=q_fac2(m6)
qq1(i7)=str_quality_1(m6)
qq2(i7)=str_quality_2(m6)
bondq4(i7)=bondq(m6)

!do i3=1,nmbond
!bondc(i7,i3)=bond_count(m6,i3)
!enddo
!do i3=1,ndet
!fvec(i7,i3)=ffvec2(m6,i3)
!enddo

if(i7.eq.strn)then

!write(9,*),'srt',(strno(m119),m119=1,i7)
!ttqlty=0
!hqlty=0
!if(nfset.gt.0)then

ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo

!if(set_number.eq.1)hqlty=ttqlty
!if(set_number.eq.1.and.nfset.eq.3)hqlty=ttqlty+noq
!print*,'hqlty=ttqlty',hqlty,ttqlty
!if(ttqlty.le.hqlty)then
!if(ttqlty1.ge.tqlty.and.ttqlty2.ge.bqlty.and.ttqlty3.ge.sqlty.or.nfset.eq.0)then
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
endif

do i4=1,100
detc(i4)=0
detc1(i4)=0
bdpstr(i4)=0
enddo

!do m119=1,i7
!do i4=1,nmbond
!detc(i4)=detc(i4)+bondc(m119,i4)
!enddo
!do i4=1,nmbond
!bdpstr(m119)=bdpstr(m119)+bondc(m119,i4)
!enddo
!enddo
!
!do i=1,nmbond
!if(detc(i).lt.6)then
!nnn=nnn+1
!endif
!enddo
!if(Ifail.eq.1.and.nnn.eq.nmbond)goto 5117
!goto 213

5117 do m119=1,i7

if(niao.eq.0)then
write(9+u1,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(bondc(m119,i4),i4=1,nmbond)
!write(9,231)(fvec(m119,i4),i4=1,ndet)
!do i4=1,ndet
!detc(fvec(m119,i4))=detc(fvec(m119,i4))+1
!enddo

!do i4=1,nmbond
!detc(i4)=detc(i4)+bondc(m119,i4)
!enddo
!do i4=1,nmbond
!bdpstr(m119)=bdpstr(m119)+bondc(m119,i4)
!enddo


enddo

!stop
!call check_ind(detc,str2,i7,nlp,nmbond,bondc)
!stop

!print*,'sourav1'
!call bond_det_vec(strno,totstr,bond_det)

call bond_det_vec(nl,i7,str2,bnd_std)

!print*,'sourav2'
!do i=1,nbd
!if(bdet(i).eq.int(bond_det))goto 778
!enddo
!nbd=nbd+1
!bdet(nbd)=int(bond_det)
!778 write(9,913)'bond_vec_det =',bond_det,Ifail,Ifailn

!if(nfset.eq.3)call Rumer_set_id(str2,i7,nl,Rid,set_num)

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif

!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

if(Rid.eq.0)write(9+u1,921)'quality value',ttqlty,'bnd_std',bnd_std
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
numb=nmbond

!write(9,232)'bond',(detc(i),i=1,numb)
!write(9,232)'det',(detc(i),i=1,tot_ndet-1)
!do i=1,tot_ndet-1
!avg=avg+detc(i)
!enddo
!avg=avg/(tot_ndet-1)
!do i=1,tot_ndet-1
!var=var+(detc(i)-avg)**2
!enddo
!var=var/(tot_ndet-1)
!std=sqrt(var)
!write(9,233)'STD_DET',std
!do i=1,1000
!detc(i)=0
!enddo
!avg=0
!var=0.0
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

!if(nfset.eq.1) write(9,*),'ovval',1.0-ovlp
!if(nfset.eq.1) goto 316
!if(set_number.eq.5000) goto 316
if(mns.eq.max_set)then
!print*,'mns',mns
!print*,'u1',u1
close(9+u1)
mns=0
u1=u1+1
!write(a,'(I1)')u1
!a=trim(a)
!print*,'a',a
outfile='structure_set_'//trim(int_num(u1))//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
goto 213
endif

i14i7=i7
m14m21=m21

do m5=m6+1,ncqs

!!!! quality 10 structures !!!!!!!
!write(*,*)'15 15 15 15 15 15'


do m20=i14i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i14i7
totstr=i14i7
m21=m14m21

totstr=totstr+1
strno(totstr)=m5
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!print*,'Ifail_test_check_ind',Ifail
call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail10',ifail,(strno(i),i=1,totstr)

if(Ifail.eq.1)then
!strno(totstr)=0
!totstr=totstr-1
goto 214
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!if(nfset.eq.1.and.totstr.gt.1) then
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
!print*,'MatLDR1'
call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 214
endif
m21=m21+1

if(m21.gt.strn) goto 316


i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m5,i3)
enddo
col(i7)=m5
!print*,'i7',i7,m5

qq(i7)=q_fac2(m5)
qq1(i7)=str_quality_1(m5)
qq2(i7)=str_quality_2(m5)
bondq4(i7)=bondq(m5)

!do i3=1,ndet
!fvec(i7,i3)=ffvec2(m5,i3)
!enddo

if(i7.eq.strn)then

!write(9,*),'srt',(strno(m119),m119=1,i7)
!ttqlty=0
!hqlty=0
!if(nfset.gt.0)then

ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo

!endif
!if(set_number.eq.1)hqlty=ttqlty
!if(set_number.eq.1.and.nfset.eq.3)hqlty=ttqlty+noq
!if(ttqlty.le.hqlty)then
!if(ttqlty1.ge.tqlty.and.ttqlty2.ge.bqlty.and.ttqlty3.ge.sqlty.or.nfset.eq.2)then
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
endif

do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(fvec(m119,i4),i4=1,ndet)
!do i4=1,ndet
!detc(fvec(m119,i4))=detc(fvec(m119,i4))+1
!enddo



enddo

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif
!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!write(9,232)'det',(detc(i),i=1,tot_ndet-1)
!do i=1,tot_ndet-1
!avg=avg+detc(i)
!enddo
!avg=avg/(tot_ndet-1)
!do i=1,tot_ndet-1
!var=var+(detc(i)-avg)**2
!enddo
!var=var/(tot_ndet-1)
!std=sqrt(var)
!write(9,233)'STD_DET',std
!do i=1,1000
!detc(i)=0
!enddo
!avg=0
!var=0.0
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!if(nfset.eq.1) write(9,*),'ovval',1.0-ovlp
!if(nfset.eq.1) goto 316
goto 214
endif
i15i7=i7
m15m21=m21

do m4=m5+1,ncqs

!!!! quality 10 structures !!!!!!!
!write(*,*)'16 16 16 16 16 16'


do m20=i15i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i15i7
totstr=i15i7
m21=m15m21

totstr=totstr+1
strno(totstr)=m4
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail10',ifail,(strno(i),i=1,totstr)

if(Ifail.eq.1)then
!strno(totstr)=0
!totstr=totstr-1
goto 215
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
!print*,'MatLDR1'
call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 215
endif
m21=m21+1

if(m21.gt.strn) goto 316


i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m4,i3)
enddo
col(i7)=m4
print*,'i7',i7,m4

qq(i7)=q_fac2(m4)
qq1(i7)=str_quality_1(m4)
qq2(i7)=str_quality_2(m4)
bondq4(i7)=bondq(m4)

!do i3=1,ndet
!fvec(i7,i3)=ffvec2(m4,i3)
!enddo

if(i7.eq.strn)then

!write(9,*),'srt',(strno(m119),m119=1,i7)
!ttqlty=0
!hqlty=0
!if(nfset.gt.0)then

ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo

!endif
!if(set_number.eq.1)hqlty=ttqlty
!if(set_number.eq.1.and.nfset.eq.3)hqlty=ttqlty+noq
!if(ttqlty.le.hqlty)then
!if(ttqlty1.ge.tqlty.and.ttqlty2.ge.bqlty.and.ttqlty3.ge.sqlty.or.nfset.eq.0)then
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
endif

do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(fvec(m119,i4),i4=1,ndet)
!do i4=1,ndet
!detc(fvec(m119,i4))=detc(fvec(m119,i4))+1
!enddo



enddo

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif
!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!write(9,232)'det',(detc(i),i=1,tot_ndet-1)
!do i=1,tot_ndet-1
!avg=avg+detc(i)
!enddo
!avg=avg/(tot_ndet-1)
!do i=1,tot_ndet-1
!var=var+(detc(i)-avg)**2
!enddo
!var=var/(tot_ndet-1)
!std=sqrt(var)
!write(9,233)'STD_DET',std
!do i=1,1000
!detc(i)=0
!enddo
!avg=0
!var=0.0
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

!if(nfset.eq.1) write(9,*),'ovval',1.0-ovlp
!if(nfset.eq.1) goto 316
goto 215
endif
i16i7=i7
m16m21=m21

!do m9=1+qstrn(10),qstrn(11)
do m3=m4+1,ncqs

!!!! quality 10 structures !!!!!!!
!write(*,*)'17 17 17 17 17 17'


do m20=i16i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i16i7
totstr=i16i7
m21=m16m21

totstr=totstr+1
strno(totstr)=m3
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail10',ifail,(strno(i),i=1,totstr)

if(Ifail.eq.1)then
!strno(totstr)=0
!totstr=totstr-1
goto 216
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
!print*,'MatLDR1'
call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 216
endif
m21=m21+1

if(m21.gt.strn) goto 316


i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m3,i3)
enddo
col(i7)=m3
!print*,'i7',i7,m3

qq(i7)=q_fac2(m3)
qq1(i7)=str_quality_1(m3)
qq2(i7)=str_quality_2(m3)
bondq4(i7)=bondq(m3)

!do i3=1,ndet
!fvec(i7,i3)=ffvec2(m3,i3)
!enddo

if(i7.eq.strn)then

!write(9,*),'srt',(strno(m119),m119=1,i7)
!ttqlty=0
!hqlty=0
!if(nfset.gt.0)then

ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo

!endif
!if(set_number.eq.1)hqlty=ttqlty
!if(set_number.eq.1.and.nfset.eq.3)hqlty=ttqlty+noq
!if(ttqlty.le.hqlty)then
!if(ttqlty1.ge.tqlty.and.ttqlty2.ge.bqlty.and.ttqlty3.ge.sqlty.or.nfset.eq.0)then
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
endif

do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(fvec(m119,i4),i4=1,ndet)
!do i4=1,ndet
!detc(fvec(m119,i4))=detc(fvec(m119,i4))+1
!enddo

enddo

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif

!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!write(9,232)'det',(detc(i),i=1,tot_ndet-1)
!do i=1,tot_ndet-1
!avg=avg+detc(i)
!enddo
!avg=avg/(tot_ndet-1)
!do i=1,tot_ndet-1
!var=var+(detc(i)-avg)**2
!enddo
!var=var/(tot_ndet-1)
!std=sqrt(var)
!write(9,233)'STD_DET',std
!do i=1,1000
!detc(i)=0
!enddo
!avg=0
!var=0.0
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!f(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!if(nfset.eq.1) write(9,*),'ovval',1.0-ovlp
!if(nfset.eq.1) goto 316
goto 216
endif
i17i7=i7
m17m21=m21

do m2=m3+1,ncqs

!!!! quality 10 structures !!!!!!!
!write(*,*)'18 18 18 18 18 18'


do m20=i17i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i17i7
totstr=i17i7
m21=m17m21

totstr=totstr+1
strno(totstr)=m2
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail10',ifail,(strno(i),i=1,totstr)

if(Ifail.eq.1)then
!strno(totstr)=0
!totstr=totstr-1
goto 217
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
!print*,'MatLDR1'
call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 217
endif
m21=m21+1

if(m21.gt.strn) goto 316


i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m2,i3)
enddo
col(i7)=m2
!print*,'i7',i7,m2

qq(i7)=q_fac2(m2)
qq1(i7)=str_quality_1(m2)
qq2(i7)=str_quality_2(m2)
bondq4(i7)=bondq(m2)

!do i3=1,ndet
!fvec(i7,i3)=ffvec2(m2,i3)
!enddo

if(i7.eq.strn)then

!write(9,*),'srt',(strno(m119),m119=1,i7)
!ttqlty=0
!hqlty=0
!if(nfset.gt.0)then

ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo

!endif
!if(set_number.eq.1)hqlty=ttqlty
!if(set_number.eq.1.and.nfset.eq.3)hqlty=ttqlty+noq
!if(ttqlty.le.hqlty)then
!if(ttqlty1.le.tqlty.and.ttqlty2.le.bqlty.and.ttqlty3.le.sqlty.or.nfset.eq.0)then
!if(ttqlty1.le.tqlty.and.ttqlty2.le.bqlty.and.ttqlty3.le.sqlty)then
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0

!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
endif

do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif

!write(9,231)(fvec(m119,i4),i4=1,ndet)
!do i4=1,ndet
!detc(fvec(m119,i4))=detc(fvec(m119,i4))+1
!enddo

enddo

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif
!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.gt.0.and.ovopt.eq.1) then
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'

if(Rid.eq.0)write(9,910)'quality value',ttqlty
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

!write(9,232)'det',(detc(i),i=1,tot_ndet-1)
!do i=1,tot_ndet-1
!avg=avg+detc(i)
!enddo
!avg=avg/(tot_ndet-1)
!do i=1,tot_ndet-1
!var=var+(detc(i)-avg)**2
!enddo
!var=var/(tot_ndet-1)
!std=sqrt(var)
!write(9,233)'STD_DET',std
!do i=1,1000
!detc(i)=0
!enddo
!avg=0
!var=0.0
!write(9,910)'qualities:',' intra_bond =',tqlty,'nn_bond =',bqlty,'sym_break=',sqlty
!write(9,911),tqlty,sqlty,bqlty
!write(9,*),'Total Quality = ',tqlty,bqlty

bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif

!if(nfset.eq.0.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp
!if(ovopt.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3) write(9,*),&
!'ovval',1.0-ovlp

if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

!if(nfset.eq.1) write(9,*),'ovlpval',1.0-ovlp
!if(nfset.eq.1) goto 316
goto 217
endif
i18i7=i7
m18m21=m21

do m1=m2+1,ncqs

!!!! quality 10 structures !!!!!!!
!write(*,*)'19 19 19 19 19 19',ifail

!print*,'i18i7,i7,m18m21,m21',i18i7,i7,m18m21,m21
if(ifail.eq.1)then
!print*,'i18i7,i7,m18m21,m21',i18i7,i7,m18m21,m21
do m20=i18i7+1,i7
strno(m20)=0
do i8=1,15
str2(m20,i8)=0
enddo
enddo

i7=i18i7
totstr=i18i7
m21=m18m21
endif

totstr=totstr+1
strno(totstr)=m1
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!print*,'Ifail_test_check_ind',Ifail
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!print*,'Ifail_test_check_ind',Ifail,Ifailn
if(Ifail.eq.1)then
goto 218
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
if(totstr.gt.1.and.ovopt.eq.vpt) then
call MatLDR('str',strno,totstr,D)

ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
if(1.0-ovlp.gt.ovval)goto 218
endif
m21=m21+1

if(m21.gt.strn) goto 316


i7=i7+1

do i3=1,nae
str2(i7,i3)=str3(m1,i3)
enddo
col(i7)=m1
!print*,'i7',i7,m10

qq(i7)=q_fac2(m1)
qq1(i7)=str_quality_1(m1)
qq2(i7)=str_quality_2(m1)
bondq4(i7)=bondq(m1)


i18i7=i7
m18m21=m21
if(i7.eq.strn)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
tqlty=tqlty+qq1(m119)
bqlty=bqlty+bondq4(m119)
sqlty=sqlty+qq2(m119)
enddo
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
set_number=set_number+1
mns=mns+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
endif



do m119=1,i7

if(niao.eq.0)then
write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
if(niao.gt.1)then
write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
if(niao.eq.1)then
write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif



enddo

call bond_det_vec(nl,i7,str2,bnd_std)

!do i=1,nbd
!if(bdet(i).eq.int(bond_det))goto 779
!enddo
!nbd=nbd+1
!bdet(nbd)=int(bond_det)
!779 write(9,913)'bond_vec_det =',bond_det,Ifail,Ifailn
!if(nfset.eq.3)call Rumer_set_id(str2,i7,nl,Rid,set_num)

if(nfset.eq.3)then
Rumwrite=1
call Rumer_set_id(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9,912)'Overlap of this set of the structures =',1.0-ovlp
endif

!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)

if(Rid.eq.0)write(9,921)'quality value',ttqlty,'bnd_std',bnd_std
if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
bqlty=0
write(9,*),'Set_number=',set_number
write(9,*),'    '
endif
if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 316
close(9+u1)
mns=0
u1=u1+1
if(u1.le.9)then
write(a,'(I1)')u1
endif
if(u1.gt.9.and.u1.lt.100)then
write(a,'(I2)')u1
endif
if(u1.gt.99.and.u1.lt.1000)then
write(a,'(I3)')u1
endif
if(u1.eq.1000)then
write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
sets. Please increase the limit'
goto 316
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

!do m20=i18i7+1,i7
!strno(m20)=0
!do i8=1,15
!str2(m20,i8)=0
!enddo
!enddo
!
!i7=i18i7
!totstr=i18i7
!m21=m18m21

goto 218
endif

!i19i7=i7
!m19m21=m21
!
!!do m9=1+qstrn(10),qstrn(11)
!do mn1=m1+1,ncqs
!
!!!!! quality 10 structures !!!!!!!
!write(*,*)'20 20 20 20 20 20'
!
!
!do m20=i19i7+1,i7
!strno(m20)=0
!do i8=1,15
!str2(m20,i8)=0
!enddo
!enddo
!
!i7=i19i7
!totstr=i19i7
!m21=m19m21
!
!
!
!!do m19=m9,ncqs
!
!
!!if(q_fac2(m19).ne.q_fac2(m9))goto 92
!totstr=totstr+1
!strno(totstr)=mn1
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!if(Ifail.eq.1)then
!goto 219
!endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!call MatLDR('str',strno,totstr,D)
!
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!if(1.0-ovlp.gt.ovval)goto 219
!endif
!m21=m21+1
!
!if(m21.gt.strn) goto 316
!
!
!i7=i7+1
!
!do i3=1,nae
!str2(i7,i3)=str3(mn1,i3)
!enddo
!col(i7)=mn1
!!print*,'i7',i7,m10
!
!qq(i7)=q_fac2(mn1)
!qq1(i7)=str_quality_1(mn1)
!qq2(i7)=str_quality_2(mn1)
!bondq4(i7)=bondq(mn1)
!
!if(i7.eq.strn)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
!enddo
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!endif
!
!do m119=1,i7
!
!if(niao.eq.0)then
!write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
!endif
!if(niao.gt.1)then
!write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
!endif
!if(niao.eq.1)then
!write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
!endif
!
!
!
!enddo
!
!if(nfset.eq.3)call Rumer_set_id(str2,i7,nl,Rid,set_num)
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!bqlty=0
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
!goto 219
!endif
!
!i20i7=i7
!m20m21=m21
!
!!do m9=1+qstrn(10),qstrn(11)
!do mn2=mn1+1,ncqs
!
!!!!! quality 10 structures !!!!!!!
!write(*,*)'21 21 21 21 21 21 21 '
!
!
!!do m20=i20i7+1,i7
!!strno(m20)=0
!!do i8=1,15
!!str2(m20,i8)=0
!!enddo
!!enddo
!!
!!i7=i20i7
!!totstr=i20i7
!!m21=m20m21
!
!
!
!!do m19=m9,ncqs
!
!
!!if(q_fac2(m19).ne.q_fac2(m9))goto 92
!totstr=totstr+1
!strno(totstr)=mn2
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!call check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!if(Ifail.eq.1)then
!goto 220
!endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!call MatLDR('str',strno,totstr,D)
!
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!if(1.0-ovlp.gt.ovval)goto 220
!endif
!m21=m21+1
!
!if(m21.gt.strn) goto 316
!
!
!i7=i7+1
!
!do i3=1,nae
!str2(i7,i3)=str3(mn2,i3)
!enddo
!col(i7)=mn2
!!print*,'i7',i7,m10
!
!qq(i7)=q_fac2(mn2)
!qq1(i7)=str_quality_1(mn2)
!qq2(i7)=str_quality_2(mn2)
!bondq4(i7)=bondq(mn2)
!
!if(i7.eq.strn)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
!enddo
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!endif
!
!do m119=1,i7
!
!if(niao.eq.0)then
!write(9,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
!endif
!if(niao.gt.1)then
!write(9,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
!endif
!if(niao.eq.1)then
!write(9,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
!endif
!
!
!
!enddo
!
!if(nfset.eq.3)call Rumer_set_id(str2,i7,nl,Rid,set_num)
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!if(Rid.eq.0)write(9,910)'quality value',ttqlty
!if(Rid.eq.1)write(9,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!bqlty=0
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!if(set_number.eq.1000)goto 316
!if(nfset.eq.0.and.ttqlty.le.ttqlty0) goto 316
!goto 220
!
!do m20=i20i7+1,i7
!strno(m20)=0
!do i8=1,15
!str2(m20,i8)=0
!enddo
!enddo
!
!i7=i20i7
!totstr=i20i7
!m21=m20m21
!
!endif
!
!
!!if(noq.eq.11)goto 210
!220 enddo
!219 enddo
218 enddo
217 enddo
216 enddo
215 enddo
214 enddo
213 enddo
212 enddo
211 enddo
210 enddo
209 enddo
208 enddo
207 enddo
206 enddo
205 enddo
204 enddo
203 enddo
202 enddo
201 enddo
200 enddo


!900 format(25I4)
!901 format(I3,x,I1,a,I1,x,25I4)
!909 format(I3,x,I1,a,I1,x,25I4)


900 format(I3,2x,I3,2x,I3,2x,I3,x,a,x,25I4)
901 format(I3,2x,I3,2x,I3,2x,I3,x,a,x,I2,a,I2,x,25I4)
909 format(I3,2x,I3,2x,I3,2x,I3,x,a,x,I3,I3,x,25I4)
910 format(a,2x,I5)
920 format(a,2x,I5,2x,a,2x,100I5)
921 format(a,2x,I5,2x,a,2x,F10.3)
911 format(15x,I3,7x,I3,7x,I3)
912 format(a,3x,F10.3)
913 format(a,3x,F10.3,2X,2I3)

231 format(50I5)
232 format(a,300I3)

close(21)

open(unit=121,file='script10',status='unknown')
write(121,111)'rm -rf','Rumer_Sets.dat'
close(121)
CALL SYSTEM ("chmod +x script10 ")
CALL SYSTEM ("./script10 ")
CALL SYSTEM ("rm script10 ")
111 format(a,x,a)

print*,'exit write_xmi_new_2',iter_counter
!316 write(9,231)(bdet(i),i=1,nbd)
316 if(nfset.eq.0)then
do i=1,i7
nstr7=nstr7+1
do j=1,nae
str5(nstr7,j)=str2(i,j)
enddo
enddo
endif
return

end subroutine write_xmi_new_2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bond_det_vec(nl,str_num,str3,bnd_std)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,i1,i2,i3,l,ll,j,nl,totbnd,totseorb,str_bnd(15000,1000),str_rad(500,20)&
,bond_sl(10000,2),num_orb(20),strno(15000),totstr,str3(2000,20),str_num,bnd_str(1000),&
avg
Double Precision::bond_vec_mat(150,150),Eig(150),bond_det
real::bnd_std,var
common /chek/totbnd,totseorb,str_bnd,str_rad,bond_sl,num_orb
common /bnd_mat/bond_vec_mat

!print*,'totstr',str_num
!do i=1,str_num
!write(9,111),'str_bnd',(str3(i,j))
!enddo
111 format (a,2x,30I3)

do i=1,totbnd
!print*,'bond_sl',(bond_sl(i,i1),i1=1,2)
ll=0
do j=1,str_num
!write(*,111),'str',(str3(j,i1),i1=1,nae)
do i2=nl*2+1,nae-nlast,2
l=0
do i1=1,2
do i3=i2,i2+1
!print*,'str3(j,i3).eq.bond_sl(i,i1)',str3(j,i3),bond_sl(i,i1)
if(str3(j,i3).eq.bond_sl(i,i1))then
l=l+1
!print*,l
goto 300
endif
enddo
300 enddo
!print*,'l',l
if(l.eq.2)then
ll=ll+1
!print*,'ll*',ll
endif
enddo
enddo
!print*,'ll',ll
bnd_str(i)=ll
enddo

avg=0
do i=1,totbnd
avg=avg+bnd_str(i)
enddo
avg=avg/totbnd

var=0
do i=1,totbnd
var=var+(bnd_str(i)-avg)**2
enddo

bnd_std=sqrt(var/totbnd)

!write(9,*)(bnd_str(i),i=1,totbnd)
!write(*,*),'bnd_str',bnd_std

!stop

!do i=1,str_num
!write(*,111),'str',(str3(i,j),j=1,nae)
!enddo
!do i=1,totbnd
!write(*,111),'bond_sl',(bond_sl(i,i1),i1=1,2)
!enddo
!
!stop

!do i=1,150
!do i1=1,150
!bond_vec_mat(i,i1)=0.0
!enddo
!enddo

!do i=1,totstr
!do i2=1,totbnd
!str_bnd(strno(i),i2)=str_bnd(strno(i),i2)+1
!enddo
!enddo

!do i=1,totstr
!do i1=i,totstr
!do i2=1,totbnd
!bond_vec_mat(i,i1)=bond_vec_mat(i,i1)+str_bnd(strno(i),i2)*str_bnd(strno(i1),i2)
!
!enddo
!enddo
!enddo
!
!
!do i=1,totstr
!do i1=i,totstr
!bond_vec_mat(i1,i)=bond_vec_mat(i,i1)
!
!enddo
!enddo
!
!!do i=1,totstr
!!write(9,111),'bond_vec_mat',(int(bond_vec_mat(i,i1)),i1=1,totstr)
!!enddo 
!
!call MatLDR('bnv',strno,totstr,Eig)
!
!!write(9,*)(Eig(i),i=1,totstr)
!
!bond_det=1.0
!do i=1,totstr
!bond_det=bond_det*Eig(i)
!enddo
!write(9,*),'bond_det',bond_det

!stop
return
end subroutine bond_det_vec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eq_dstr_set(tns,lnp,npstr,str2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,&
l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,j
integer::tns,lnp,npstr,fail,strsln(1000),str2(15000,20),num_str,numpstr
common/nst/num_str,numpstr

numpstr=npstr
num_str=0

print*,'enter eq_dstr_set'
!write(9,*),'tns,lnp,npstr',tns,lnp,npstr
!call check_str_bond(lnp,tns,str2,npstr)

do i=1,tns
write(*,200)i,(str2(i,i1),i1=1,nae)
enddo
200 format (30I5)

!stop
do i1=1,tns
do i=1,200
strsln(i)=0
enddo
j=0
j=j+1
strsln(j)=i1

l1=j
do i2=i1+1,tns
do i=l1+1,j
strsln(i)=0
enddo
j=l1
j=j+1
strsln(j)=i2
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 102

l2=j

do i3=i2+1,tns
do i=l2+1,j
strsln(i)=0
enddo
j=l2
j=j+1
strsln(j)=i3
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 103

l3=j


do i4=i3+1,tns
do i=l3+1,j
strsln(i)=0
enddo
j=l3
j=j+1
strsln(j)=i4
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 104

l4=j


do i5=i4+1,tns
do i=l4+1,j
strsln(i)=0
enddo
j=l4
j=j+1
strsln(j)=i5 
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 105

l5=j


do i6=i5+1,tns
do i=l5+1,j
strsln(i)=0
enddo
j=l5
j=j+1
strsln(j)=i6
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 106

l6=j


do i7=i6+1,tns
do i=l6+1,j
strsln(i)=0
enddo
j=l6
j=j+1
strsln(j)=i7
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 107

l7=j

do i8=i7+1,tns
do i=l7+1,j
strsln(i)=0
enddo
j=l7
j=j+1
strsln(j)=i8
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 108

l8=j

do i9=i8+1,tns
do i=l8+1,j
strsln(i)=0
enddo
j=l8
j=j+1
strsln(j)=i9
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 109

l9=j

do i10=i9+1,tns
do i=l9+1,j
strsln(i)=0
enddo
j=l9
j=j+1
strsln(j)=i10
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 110

l10=j

do i11=i10+1,tns
do i=l10+1,j
strsln(i)=0
enddo
j=l10
j=j+1
strsln(j)=i11
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 111

l11=j

do i12=i11+1,tns
do i=l11+1,j
strsln(i)=0
enddo
j=l11
j=j+1
strsln(j)=i12
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 112

l12=j

do i13=i12+1,tns
do i=l12+1,j
strsln(i)=0
enddo
j=l12
j=j+1
strsln(j)=i13
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 113

l13=j

do i14=i13+1,tns
do i=l13+1,j
strsln(i)=0
enddo
j=l13
j=j+1
strsln(j)=i14
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 114

l14=j

do i15=i14+1,tns
do i=l14+1,j
strsln(i)=0
enddo
j=l14
j=j+1
strsln(j)=i15
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 115

l15=j

do i16=i15+1,tns
do i=l15+1,j
strsln(i)=0
enddo
j=l15
j=j+1
strsln(j)=i16
call eq_dst_check(j,strsln,fail,lnp,str2)
if (fail.eq.1) goto 116

l16=j

116 enddo
115 enddo
114 enddo
113 enddo
112 enddo
111 enddo
110 enddo
109 enddo
108 enddo
107 enddo
106 enddo
105 enddo
104 enddo
103 enddo
102 enddo
enddo


print*,'exit eq_dstr_set'
return
end subroutine eq_dstr_set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eq_dst_check(nts,strsln,fail,nl,str2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
integer::totbnd,totseorb,str_bnd(15000,1000),str_rad(500,20),bond_sl(10000,2),num_orb(20)
integer::i,i1,j,m,m1,l,add_edst(500),n1,n2,n3,quality_fac(15000),str_quality_1(15000),&
str_quality_2(15000),nssym,qulsym(15000),symq(15000),tqlty,bqlty,sqlty,tnqs,sigsym(15000)&
,tnqs_sig,fail,strsln(1000),nl,tnbd,nts,nbd,str2(15000,20),bondq(15000),strflg,num_str,numpstr,det_inv
integer::str3(15000,20),fvec(15000,1000),strsln1(1000)
real*8::factorial
common /chek/totbnd,totseorb,str_bnd,str_rad,bond_sl,num_orb
common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig
common/nst/num_str,numpstr


!write(9,*),'**********************************'
!do i=1,nts
!write(*,231),'str',(str2(i,j),j=1,nae)
!enddo
!print*,'totbnd',totbnd
!print*,'enter eq_dstr_set'
!if(strflg.eq.0)num_str=0
!print*,'strsln',(strsln(i),i=1,nts)
!do i=1,nts
!write(9,231),'str_bnd',(str_bnd(strsln(i),j),j=1,totbnd)
!enddo

!do i=1,500
!add_edst(i)=0
!enddo



fail=0
do l=1,totbnd
add_edst(l)=0
do i=1,nts
!write(9,*),'i',i,str_bnd(strsln(i),l)
add_edst(l)=add_edst(l)+str_bnd(strsln(i),l)
enddo
!write(9,*),'add_edst',l,add_edst(l)
enddo
!write(9,231),'strsln',(strsln(i),i=1,nts)
!write(9,231),'add_edst',(add_edst(i),i=1,totbnd)

nbd=((nao-nl-nlast)/2)*numpstr
n1=nao-nl
n2=nao-nl-2
!print*,'n1,n2',n1,n2
tnbd=(int(factorial(n1)/(factorial(n2)))/2)
!print*,'tnbd',tnbd,nbd
if(nbd.gt.tnbd)n3=nbd/tnbd
if(nbd.le.tnbd)n3=tnbd/nbd
do i=1,totbnd
if(add_edst(i).gt.n3)fail=1
enddo

!print*,'nbd,n1,n2,tnbd,n3,nao',nbd,n1,n2,tnbd,n3,nao,fail,j,numpstr
if(nts.eq.numpstr.and.fail.eq.0)then
do i=1,numpstr
strsln1(i)=i
do i1=1,nae
str3(i,i1)=str2(strsln(i),i1)
enddo
enddo

!do i=1,numpstr
!write(9,102)(str3(i,i1),i1=1,nae)
!enddo
102 format (50I5)

!write(9,*)'calling vector_rep'
call vector_rep(nl,str3,numpstr,fvec)
!write(9,*)'ending vector_rep'
!write(9,*),'tnbd,nbd',tnbd,nbd
if(tnbd.ne.nbd)call mat_ind(nl,nts,numpstr,strsln1,fail,det_inv)
!call mat_ind(nl,nts,numpstr,strsln,fail,det_inv)
if(fail.eq.0)then
!write(9,*),'nbd,n1,n2,tnbd,n3,nao',nbd,n1,n2,tnbd,n3,nao,fail,j,numpstr
!do i=1,nts
!write(9,231),'str',(str2(i,j),j=1,nae)
!enddo
num_str=num_str+1
!print*,'j,numpstr,num_str',j,numpstr,num_str
!write(9,231),'j,numpstr,num_str',j,numpstr,num_str
do m=1,numpstr
if(niao.eq.0)then
write(9,900),str_quality_1(strsln(m)),bondq(strsln(m)),str_quality_2(strsln(m)),'|',(str2(strsln(m),m1),m1=1,nae)
endif
if(niao.gt.1)then
write(9,901),str_quality_1(strsln(m)),bondq(strsln(m)),str_quality_2(strsln(m)),'|',1,':',niao,(str2(strsln(m),m1),m1=1,nae)
endif
if(niao.eq.1)then
write(9,909),str_quality_1(strsln(m)),bondq(strsln(m)),str_quality_2(strsln(m)),'|',1,1,(str2(strsln(m),m1),m1=1,nae)
endif
enddo

fail=1
write(9,*),''
write(9,*),'set number',num_str

write(9,*),''
write(9,*),''
endif
endif
if(num_str.eq.1000)stop
900 format(I3,2x,I3,2x,I3,2x,a,x,25I4)
901 format(I3,2x,I3,2x,I3,2x,a,x,I2,a,I2,x,25I4)
909 format(I3,2x,I3,2x,I3,2x,a,x,I3,I3,x,25I4)
231 format(a,2x,50I3)
print*,'exit eq_dstr_set'
return
end subroutine eq_dst_check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_ind(str3,totstr,nl,ncqs,strno,Ifailn,ffvec2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::j,i,ii,iii,str3(15000,20),d,e,elporb,str_rad_1(100),str_rad_2(100),ffvec2(15000,1000),&
Ifailn,i7,wig2,na,i1,i2,i3,totstr,ncqs,strno(15000),str_bnd_1(1000),str_bnd_2(1000),sing(1000,100),&
str_sl(1000),new_mult,nl,totbnd,totseorb,str_bnd(15000,1000),str_rad(500,20),bond_sl(10000,2),&
npstr,num_orb(20),val(20),val1(20),i4,i5,i6,l,l1,rul,bdsl(20),perm(20),Ifailn1
real*8::factorial
character(9)::sbrtn

common /chek/totbnd,totseorb,str_bnd,str_rad,bond_sl,num_orb
!character(5)::a

Ifailn=0

!do i=1,totbnd
!print*,(bond_sl(i,i1),i1=1,2)
!enddo

e=nae-nl*2
call wigner(e,wig2)
npstr=wig2

!print*,'totbnd',totbnd
iii=0

!do i=1,totstr
!write(*,100),'str3',(str3(strno(i),ii),ii=1,nae),strno(i)
!!write(*,100),'fvec',(ffvec2(strno(i),ii)*sing(i,ii),ii=1,32)
!enddo

if(nlast.ne.0)then
do ii=1,totseorb
do i=1,totstr
str_rad_1(ii)=str_rad_1(ii)+str_rad(strno(i),ii)
enddo
enddo
d=0
do ii=1,totseorb
if(str_rad_1(ii).gt.d)d=str_rad_1(ii)
enddo
if(d.lt.3) goto 121
new_mult=mult
e=nae-1-nl*2
mult=mult-1
call wigner(e,wig2)
mult=new_mult
do ii=1,totseorb
if(str_rad_1(ii).gt.wig2)then
Ifailn=1
goto 121
endif
enddo

do ii=3,d-1
do i=1,totseorb
if(ii.eq.str_rad_1(i))then
do i1=1,totstr
if(str_rad(strno(i1),i).eq.1)then
iii=iii+1
str_sl(iii)=strno(i1)
endif
enddo

do i2=1,totseorb
do i1=1,iii
str_rad_2(i2)=str_rad_2(i2)+str_rad(str_sl(i1),i2)
enddo
enddo

i3=0
do i2=1,totseorb
if(str_rad_2(i2).eq.ii)then
i3=i3+1
endif
enddo

new_mult=mult
e=nae-nl*2-i3
mult=(nlast-i3)+1
call wigner(e,wig2)
mult=new_mult
if(ii.gt.wig2)then
Ifailn=1
goto 121
endif
endif
enddo
enddo


endif



do ii=1,totbnd
!print*,'ii',ii
str_bnd_1(ii)=0
do i=1,totstr
str_bnd_1(ii)=str_bnd_1(ii)+str_bnd(strno(i),ii)
enddo
!print*,'str_bnd',str_bnd_1(ii),ii
enddo

d=0
do ii=1,totbnd
if(str_bnd_1(ii).gt.d)d=str_bnd_1(ii)
enddo

!if(totstr.lt.0)then
!do i1=1,nae-nl*2
!do i2=i1+1,nae-nl*2
!do i3=i2+1,nae-nl*2
!val(1)=num_orb(i1)
!val(2)=num_orb(i2)
!val(3)=num_orb(i3)
!
!print*,'orbs',(val(i6),i6=1,3)
!l1=0
!do i4=1,totbnd
!l=0
!do i5=1,2
!do i6=1,3
!if(bond_sl(i4,i5).eq.val(i6))then
!l=l+1
!endif
!enddo
!if(l.eq.2)then
!l1=l1+1
!val1(l1)=i4
!endif
!enddo
!enddo
!
!!print*,'orbnum',(val1(i6),i6=1,l1)
!l=0
!do i4=1,l1
!l=l+str_bnd_1(val1(i4))
!enddo
!rul=int(totstr/3)
!!print*,'lllll',l
!if(l.gt.rul*2)then
!Ifailn=1
!
!print*,'orbnum',rul,l,'>',(val1(i6),i6=1,l1)
!goto 121
!endif
!enddo
!enddo
!enddo
!
!!stop
!
!endif
!!l1=npstr/3



do ii=1,totstr
!write(*,100),'str',(str3(strno(ii),i),i=1,nae)
write(*,100),'str_bnd',strno(ii),(str_bnd(strno(ii),i),i=1,totbnd)
enddo
write(*,100),'add_val',(str_bnd_1(i),i=1,totbnd)
print*,'ddd',d,totbnd
100 format(a,50I3)

if(d.lt.3) goto 121
e=nae-2-nl*2
call wigner(e,wig2)
!print*,'wig2',wig2
do ii=1,totbnd
if(str_bnd_1(ii).gt.wig2)then
Ifailn=1
goto 121
endif
enddo



do i=1,totbnd
if(str_bnd_1(i).gt.2.and.str_bnd(strno(totstr),i).eq.1)then
iii=0
do i1=1,totstr
if(str_bnd(strno(i1),i).eq.1)then
iii=iii+1
str_sl(iii)=strno(i1)
print*,'str_sl',i,'>',iii,str_sl(iii)
endif
enddo

call all_set_gen(str_sl,iii,nl,i,Ifailn)
if(Ifailn.eq.1)goto 121
endif
enddo

!stop

ii=0
do i=1,totbnd
if(str_bnd(strno(totstr),i).eq.1)then
ii=ii+1
bdsl(ii)=i
endif
enddo

sbrtn='check_ind'
do j=2,int(ii/2)

do i1=1,ii
perm(1)=bdsl(i1)

do i2=i1+1,ii
perm(2)=bdsl(i2)
if(j.eq.2)then
print*,(perm(i),i=1,j),j
!stop
call check3(j,perm,strno,totstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 122
goto 101
endif

do i3=i2+1,ii
perm(3)=bdsl(i3)
if(j.eq.3)then
call check3(j,perm,strno,totstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 122
goto 102
endif

do i4=i3+1,ii
perm(4)=bdsl(i4)
if(j.eq.4)then
call check3(j,perm,strno,totstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 122
goto 103
endif

do i5=i4+1,ii
perm(5)=bdsl(i5)
if(j.eq.5)then
call check3(j,perm,strno,totstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 122
goto 104
endif

do i6=i5+1,ii
perm(6)=bdsl(i6)
if(j.eq.6)then
call check3(j,perm,strno,totstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 122
goto 105
endif

do i7=i6+1,ii
perm(7)=bdsl(i7)
if(j.eq.7)then
call check3(j,perm,strno,totstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 122
goto 106
endif

106 enddo
105 enddo
104 enddo
103 enddo
102 enddo
101 enddo
    enddo
 enddo

122 Ifailn=Ifailn1



121 print*,'Ifailn',Ifailn
print*,'ssssssssssssssssssssssssssssssssssssssssssssssssssssss'
!if(Ifailn.eq.1)stop
return
end subroutine check_ind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check3(e,perm,str_sl,totstr,nl,ifailn1,sbrtn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,e,a,nnstr,ifailn1,perm(20),wig2,nl,str_sl(1000),&
totbnd,totseorb,str_bnd(15000,1000),str_rad(500,20),bond_sl(10000,2),wig3,e2,e1,bnd(20),new_bd(100,2)
integer::totstr,bd_sl(20),num_orb(20)
character(9)::sbrtn

common /chek/totbnd,totseorb,str_bnd,str_rad,bond_sl,num_orb

print*,'***',sbrtn,'***'
do i=1,totstr
write(*,115),'check3',(str_bnd(str_sl(i),i1),i1=1,totbnd)
enddo
print*,'perm',(perm(i),i=1,e)
print*,'sl_num',(str_sl(i),i=1,totstr)
115 format(a,x,50I3)
ifailn1=0
i2=0
do i=1,e
do i1=1,2
i2=i2+1
bnd(i2)=bond_sl(perm(i),i1)
enddo
enddo
print*,'bnd',(bnd(i),i=1,i2)

a=0
do i=1,i2
!new_bd(1)=bnd(i)
do i1=i+1,i2
a=a+1
new_bd(a,1)=bnd(i)
new_bd(a,2)=bnd(i1)
enddo
enddo

do i=1,a
print*,'new_bd',(new_bd(i,i1),i1=1,2)
enddo

do i=1,a
do i3=1,totbnd
i6=0
do i4=1,2
do i5=1,2
if(new_bd(i,i4).eq.bond_sl(i3,i5))then
i6=i6+1
endif
enddo
enddo
if(i6.eq.2)then
bd_sl(i)=i3
goto 120
endif
enddo
120 enddo
print*,'bd_sl',(bd_sl(i),i=1,a)


i8=0
do i=1,totstr
i7=0
do i5=1,a
if(str_bnd(str_sl(i),bd_sl(i5)).eq.1)then
print*,'bd_sl(i5)',bd_sl(i5),i,i5
i7=i7+1
endif
enddo
if(i7.eq.e)then
print*,'i',i
i8=i8+1
endif
enddo

if(sbrtn.eq.'all_set'.or.sbrtn.eq.'check_2')then
e1=nae-nl*2-e*2-2
else
e1=nae-nl*2-e*2
endif
e2=e*2
print*,'e1,e2',e1,e2
call wigner(e1,wig2)
call wigner(e2,wig3)
if(i8.gt.wig3*wig2)ifailn1=1
print*,'ifailn1',ifailn1,nl,nae,e,wig2,wig3,i8
!stop

print*,'exit check3'
return
end subroutine check3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_ind_2(permut,nnstr,nl,n,Ifailn)

use commondat
implicit none

integer::i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,e,a,nnstr,Ifailn,permut(100),num_orb(20),wig2,nl
integer::str_bnd_2(1000),totbnd,totseorb,str_bnd(15000,1000),str_rad(500,20),bond_sl(10000,2)
integer::bdsl(20),perm(20),n,ii,e1,e2,wig3,j,Ifailn1
character(9)::sbrtn

common /chek/totbnd,totseorb,str_bnd,str_rad,bond_sl,num_orb

print*,'enter check_ind_2'
print*,'permut',(permut(i),i=1,nnstr)
Ifailn=0

do i2=1,totbnd
str_bnd_2(i2)=0
do i1=1,nnstr
str_bnd_2(i2)=str_bnd_2(i2)+str_bnd(permut(i1),i2)
enddo
enddo
!
!do i2=1,iii
!!write(*,100),'str',(str3(strno(ii),i),i=1,nae)
!write(*,100),'str_bnd_2',(str_bnd(str_sl(i2),i1),i1=1,totbnd)
!enddo
!write(*,100),'add_val_2',(str_bnd_2(i1),i1=1,totbnd)
!
i3=0
do i2=1,totbnd
if(str_bnd_2(i2).eq.nnstr)then
i3=i3+1
endif
enddo

e=nae-nl*2-i3*2
call wigner(e,wig2)
!print*,'wig2ll',wig2,i3
if(nnstr.gt.wig2)then
Ifailn=1
endif
print*,'Ifailn',Ifailn
if(Ifailn.eq.1)goto 100

!ii=0
!do i=1,totbnd
!if(i.eq.n)goto 400
!if(str_bnd(permut(nnstr),i).eq.1)then
!ii=ii+1
!bdsl(ii)=i
!endif
!400 enddo
!print*,(bdsl(i),i=1,ii)
!!stop
!sbrtn='check_2'
!!if(ii)
!do j=2,int(ii/2)
!e1=nae-nl*2-j*2-2
!e2=j*2
!print*,'e1,e2',e1,e2
!call wigner2(e1,wig2)
!call wigner2(e2,wig3)
!
!print*,wig2*wig3,nnstr
!if(wig2*wig3.gt.nnstr)goto 300
!
!do i1=1,ii
!perm(1)=bdsl(i1)
!
!do i2=i1+1,ii
!perm(2)=bdsl(i2)
!if(j.eq.2)then
!!print*,(perm(i),i=1,j),j
!!stop
!call check3(j,perm,permut,nnstr,nl,Ifailn1,sbrtn)
!if(Ifailn1.eq.1)goto 122
!goto 201
!endif
!
!do i3=i2+1,ii
!perm(3)=bdsl(i3)
!if(j.eq.3)then
!call check3(j,perm,permut,nnstr,nl,Ifailn1,sbrtn)
!if(Ifailn1.eq.1)goto 122
!goto 202
!endif
!
!do i4=i3+1,ii
!perm(4)=bdsl(i4)
!if(j.eq.4)then
!call check3(j,perm,permut,nnstr,nl,Ifailn1,sbrtn)
!if(Ifailn1.eq.1)goto 122
!goto 203
!endif
!
!do i5=i4+1,ii
!perm(5)=bdsl(i5)
!if(j.eq.5)then
!call check3(j,perm,permut,nnstr,nl,Ifailn1,sbrtn)
!if(Ifailn1.eq.1)goto 122
!goto 204
!endif
!
!do i6=i5+1,ii
!perm(6)=bdsl(i6)
!if(j.eq.6)then
!call check3(j,perm,permut,nnstr,nl,Ifailn1,sbrtn)
!if(Ifailn1.eq.1)goto 122
!goto 205
!endif
!
!do i7=i6+1,ii
!perm(7)=bdsl(i7)
!if(j.eq.7)then
!call check3(j,perm,permut,nnstr,nl,Ifailn1,sbrtn)
!if(Ifailn1.eq.1)goto 122
!goto 206
!endif
!
!206 enddo
!205 enddo
!204 enddo
!203 enddo
!202 enddo
!201 enddo
!    enddo
!300 enddo
!
!122 Ifailn=Ifailn1

!print*,'Ifailn',Ifailn,wig2,i3
100 print*,'exit check_ind_2'
return
end subroutine check_ind_2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine all_set_gen(str_sl,nnstr,nl,n,Ifailn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,ii,elporb,nl,a,nnstr,Ifailn,permut(100),&
str_sl(1000),wig2,l,k,bdsl(20),perm(20),j,n
integer::e1,e2,wig3,Ifailn1
character(9)::sbrtn
integer::totbnd,totseorb,str_bnd(15000,1000),str_rad(500,20),bond_sl(10000,2),num_orb(20)

common /chek/totbnd,totseorb,str_bnd,str_rad,bond_sl,num_orb
print*,'enter all_set_gen'
print*,'str_sl',(str_sl(i),i=1,nnstr)



ii=0
do i=1,totbnd
if(i.eq.n)goto 400
if(str_bnd(str_sl(nnstr),i).eq.1)then
ii=ii+1
bdsl(ii)=i
endif
400 enddo
print*,(bdsl(i),i=1,ii)
!stop
sbrtn='all_set'
!if(ii)
do j=2,int(ii/2)
e1=nae-nl*2-j*2-2
e2=j*2
print*,'e1,e2',e1,e2
call wigner(e1,wig2)
call wigner(e2,wig3)

print*,wig2*wig3,nnstr
if(wig2*wig3.gt.nnstr)goto 300

do i1=1,ii
perm(1)=bdsl(i1)

do i2=i1+1,ii
perm(2)=bdsl(i2)
if(j.eq.2)then
!print*,(perm(i),i=1,j),j
!stop
call check3(j,perm,str_sl,nnstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 500
goto 201
endif

do i3=i2+1,ii
perm(3)=bdsl(i3)
if(j.eq.3)then
call check3(j,perm,str_sl,nnstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 500
goto 202
endif

do i4=i3+1,ii
perm(4)=bdsl(i4)
if(j.eq.4)then
call check3(j,perm,str_sl,nnstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 500
goto 203
endif

do i5=i4+1,ii
perm(5)=bdsl(i5)
if(j.eq.5)then
call check3(j,perm,str_sl,nnstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 500
goto 204
endif

do i6=i5+1,ii
perm(6)=bdsl(i6)
if(j.eq.6)then
call check3(j,perm,str_sl,nnstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 500
goto 205
endif

do i7=i6+1,ii
perm(7)=bdsl(i7)
if(j.eq.7)then
call check3(j,perm,str_sl,nnstr,nl,Ifailn1,sbrtn)
if(Ifailn1.eq.1)goto 500
goto 206
endif

206 enddo
205 enddo
204 enddo
203 enddo
202 enddo
201 enddo
    enddo
300 enddo









a=((nae-nl*2)-mod(nae-nl*2,2))/2
!print*,'aa',a,nae,mod(nae-nl*2,2)
    
do i=1,a-2
elporb=nae-nl*2-i*2
!print*,'elporb',elporb
call wigner(elporb,wig2)
!print*,'wig2',wig2,elporb,nnstr
if(wig2.lt.nnstr+1)then

!print*,'sourav'
do i1=1,nnstr-1
permut(1)=str_sl(i1)
if(wig2.eq.1)then
permut(2)=str_sl(nnstr)
l=2
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 100
endif

do i2=i1+1,nnstr-1
permut(2)=str_sl(i2)
if(wig2.eq.2)then
permut(3)=str_sl(nnstr)
l=3
!print*,'permut',(permut(k),k=1,3)
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 101
endif


do i3=i2+1,nnstr-1
permut(3)=str_sl(i3)
if(wig2.eq.3)then
permut(4)=str_sl(nnstr)
l=4
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 102
endif

do i4=i3+1,nnstr-1
permut(4)=str_sl(i4)
if(wig2.eq.4)then
permut(5)=str_sl(nnstr)
l=5
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 103
endif

do i5=i4+1,nnstr-1
permut(5)=str_sl(i5)
if(wig2.eq.5)then
permut(6)=str_sl(nnstr)
l=6
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 104
endif

do i6=i5+1,nnstr-1
permut(6)=str_sl(i6)
if(wig2.eq.6)then
permut(7)=str_sl(nnstr)
l=7
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 105
endif

do i7=i6+1,nnstr-1
permut(7)=str_sl(i7)
if(wig2.eq.7)then
permut(8)=str_sl(nnstr)
l=8
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 106
endif

do i8=i7+1,nnstr-1
permut(8)=str_sl(i8)
if(wig2.eq.8)then
permut(9)=str_sl(nnstr)
l=9
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 107
endif

do i9=i8+1,nnstr-1
permut(9)=str_sl(i9)
if(wig2.eq.9)then
permut(10)=str_sl(nnstr)
l=10
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 108
endif

do i10=i9+1,nnstr-1
permut(10)=str_sl(i10)
if(wig2.eq.10)then
permut(11)=str_sl(nnstr)
l=11
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 109
endif

do i11=i10+1,nnstr-1
permut(11)=str_sl(i11)
if(wig2.eq.11)then
permut(12)=str_sl(nnstr)
l=12
call check_ind_2(permut,l,nl,n,Ifailn)
if(Ifailn.eq.1)goto 600
goto 110
endif

110 enddo
109 enddo
108 enddo
107 enddo
106 enddo
105 enddo
104 enddo
103 enddo
102 enddo
101 enddo
100 enddo


endif
enddo
!500 stop
500 Ifailn=Ifailn1
print*,'exit all_set_gen'
600 return
end subroutine all_set_gen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Rumer_set_id(str2,i7,lonep,Rid,set_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/qfacRumid/qq1,bondq4,qq2,qq,rumset

integer::i,ii,iii,iiii,i1,i2,i3,i4,i5,i6,i7,j,Rid,orb1,orb2,sc,orbs(500),str2(2000,20),lonep,rumstrset(500,20),set_num(100)
integer::qq1(5000),qq2(5000),qq(5000),bondq4(5000),m119,m20,rumset
character(5)::a


!print*,'totrum',totrum,lonep,i7
j=0
!do i=1,i7
!print*,(str2(i,i1),i1=1,nae)
!enddo
!write(23,*)nae,lonep,(perm(i),i=1,nae-lonep*2)
!write(*,*)nae,lonep,(perm(i),i=1,nae-lonep*2)
!stop

Rid=0
rewind(31)
do i=1,totrum
iiii=0
read(31,*)
 do i1=1,i7
read(31,*)(rumstrset(i1,i2),i2=1,nae)
!print*,'RUM',(rumstrset(i1,i2),i2=1,nae)
 enddo
  do i1=1,i7
   do i2=1,i7
iii=0
    do i3=lonep*2+1,nae-nlast,2
     do i4=lonep*2+1,nae-nlast,2
ii=0
      do i5=i3,i3+1
      do i6=i4,i4+1
!print*,str2(i1,i6),rumstrset(i2,i5)
if(rumstrset(i2,i5).eq.str2(i1,i6))then
ii=ii+1
goto 249
endif
      enddo
249   enddo
if(ii.eq.2)then
iii=iii+ii
goto 250
endif
     enddo
250 enddo
if(iii.eq.nae-lonep*2-nlast)then
iiii=iiii+1
goto 251
endif
   enddo
251 enddo
if(iiii.eq.i7)then
j=j+1
Rid=1
set_num(j)=i
endif
252 enddo

nrs=j

if(Rid.eq.1.and.Rumwrite.eq.1)then
rumset=rumset+1
!write(*,*)'set number',rumset,(perm(i),i=1,nae-lonep*2)
!stop
do m119=1,i7

if(niao.eq.0)then
write(23,900),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',(str2(m119,m20),m20=1,nae)
endif
if(niao.gt.1)then
write(23,901),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,':',niao,(str2(m119,m20),m20=1,nae)
endif
if(niao.eq.1)then
write(23,909),qq1(m119),bondq4(m119),qq2(m119),qq(m119),'|',1,1,(str2(m119,m20),m20=1,nae)
endif
enddo
write(23,*),'Set_number=',rumset
write(23,*)
endif


900 format(I3,2x,I3,2x,I3,2x,I3,x,a,x,25I4)
901 format(I3,2x,I3,2x,I3,2x,I3,x,a,x,I2,a,I2,x,25I4)
909 format(I3,2x,I3,2x,I3,2x,I3,x,a,x,I3,I3,x,25I4)

!300 rewind(31)
300 return
end subroutine Rumer_set_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine charge_cntr(str2,nl,strnum)
use commondat
implicit none

integer::i,i1,i2,i3,i4,iii,j,j1,nl,str2(15000,20),strnum,nn(10)
real*8::cent_coord(5000,3),coordx(100),coordy(100),coordz(100),chrg_cntrx,chrg_cntry,chrg_cntrz,charge_p_cntr(50),&
chrg_cntr(1000,3)

common/coordinate/coordx,coordy,coordz

!i4=0

do i=1,atom
print*,'coord',coordx(i),coordy(i),coordz(i)
enddo
do i=1,strnum
j1=0
print*,'str2',(str2(i,i1),i1=1,nae)
do i1=1+nl*2,(nae-nlast),2
iii=0
nn(1)=0
nn(2)=0
do i2=i1,i1+1
do i3=1,atom
do i4=1,atn(active_atoms(i3))
if(str2(i,i2).eq.atoset(active_atoms(i3),i4))then
iii=iii+1
if(mod(i2,2).eq.0)nn(2)=active_atoms(i3)
if(mod(i2,2).eq.1)nn(1)=active_atoms(i3)
goto 101
endif
enddo
enddo
101 enddo
!if(nn(2).eq.nn(1))goto 202
if(iii.ne.2) goto 202
print*,'nn(1),nn(2)',i,nn(1),nn(2)

j1=j1+1
cent_coord(j1,1)=(coordx(nn(1))+coordx(nn(2)))/2.0
cent_coord(j1,2)=(coordy(nn(1))+coordy(nn(2)))/2.0
cent_coord(j1,3)=(coordz(nn(1))+coordz(nn(2)))/2.0
charge_p_cntr(j1)=2.0

202 enddo
!if (nlast.ne.0)then
!do i1=nae-nlast+1,nae
!do i3=1,atom
!do i4=1,atn(active_atoms(i3))
!if(str2(i,i1).eq.atoset(active_atoms(i3),i4))then
!nn(1)=active_atoms(i3)
!print*,'radical_atom',nn(1)
!j1=j1+1
!cent_coord(j1,1)=coordx(nn(1))
!cent_coord(j1,2)=coordy(nn(1))
!cent_coord(j1,3)=coordz(nn(1))
!charge_p_cntr(j1)=1.0
!endif
!enddo
!enddo
!enddo
!endif

do i1=1,j1
print*,'cent_coord',(cent_coord(i1,i2),i2=1,3)
enddo
chrg_cntrx=0
chrg_cntry=0
chrg_cntrz=0
do i1=1,j1
chrg_cntrx=chrg_cntrx+charge_p_cntr(i1)*cent_coord(i1,1)
chrg_cntry=chrg_cntry+charge_p_cntr(i1)*cent_coord(i1,2)
chrg_cntrz=chrg_cntrz+charge_p_cntr(i1)*cent_coord(i1,3)
enddo
chrg_cntr(i,1)=chrg_cntrx/(nae-nl*2)
chrg_cntr(i,2)=chrg_cntry/(nae-nl*2)
chrg_cntr(i,3)=chrg_cntrz/(nae-nl*2)
print*,'chrg_cntr',(chrg_cntr(i,j),j=1,3)
enddo


!stop

return
end subroutine charge_cntr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine str_selection(astr,nl,m4,perm_nstr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig

integer::i,j,i1,i2,i3,i4,i5,i6,i7,ii,nl,str1_cnt,cnt,m4,m8,m9,m10,lp_cnt,tnqs,&
wig2,elporb,perm_nstr,total_str,mbondq(15000),q_fac(15000),q_fac1(15000),ijk,n,nssym,rep
integer::str_cnt(15000),astr(15000,20),str1(15000,20),str2(15000,20),lps(50),sigsym(15000),tnqs_sig&
,quality_fac(15000),str_quality_1(15000),str_quality_2(15000),strno(1000)&
,rumer(15000),rumer_rad(15000),bondq(15000),fvec(15000,1000),qulsymm(15000),&
symq(15000),symqq(15000),&
tndet,tqlty,tqlty1,bqlty,sqlty,bqlty1,sqlty1,qulsym(15000)
real*8::factorial,dfactorial,symsc(15000)
character(5)::rumstr


print*,'enter str_selection'
print*,'ovopt',ovopt,nfset
tqlty1=0
bqlty1=0
sqlty1=0
!do i1=1,m4
!write(*,305),i1,(astr(i1,i2),i2=1,nae)
!enddo
!print*,perm_nstr,nl,m4

!**********************************************************************************************************
!!!!!!!!!!!!!!!! If number of permisible structures and available structures are same START !!!!!!!!!!!!!!!
!**********************************************************************************************************

if(perm_nstr.eq.m4)then
call rumer_structures(nl,astr,m4,rumer,rumer_rad)
call quality_factor(nl,astr,m4,quality_fac,str_quality_1,str_quality_2,bondq)
!call intra_bond_factor(nl,astr,m4,str_quality_1)
!call symm_break_factor(nl,astr,m4,str_quality_2)
!call symmetry_cal(nl,str1,str1_cnt,symsc,symq,nssym)


write(7,*),'all',m4,' structures are permisible among below '
!n=0
!do ijk=nssym,1,-1
do m8=1,m4
!if(symq(m8).eq.ijk)then
!n=n+1
if(rumer(m8)*rumer_rad(m8).eq.1)rumstr='R'
if(rumer(m8)*rumer_rad(m8).eq.0)rumstr='-'
write(7,301),'str',m8,')','[',str_quality_1(m8),',',str_quality_2(m8),']','{',quality_fac(m8),'}',&
rumstr,(astr(m8,m9),m9=1,nae)
!write(7,305)(fvec(m8,m9),m9=1,ndet)
!write(*,305)(str1(m8,m9),m9=1,nae)
!endif
enddo
!enddo

301 format(a,x,I5,a,x,a,I3,a,I3,x,a,a,I3,x,a,x,a,x,30I5)
305 format(30I5)

!if(nnnatom.ne.0)then
if(input_flg.eq.1)call nnat_bond_cal(nl,astr,m4,bondq)
if(input_flg.eq.0)call nnat_bond_cal_2(nl,astr,m4,bondq)
!endif
tqlty=0
bqlty=0
sqlty=0
do m8=1,m4
if(niao.eq.0)then
write(9,900),str_quality_1(m8),bondq(m8),str_quality_2(m8),'|',(astr(m8,m9),m9=1,nae)
endif
!if(niao.gt.1.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9,901),str_quality_1(m8),bondq(m8),str_quality_2(m8),'|',1,':',niao,(astr(m8,m9),m9=1,nae)
endif
!if(niao.eq.1.and.nnnatom.ne.0)then
if(niao.eq.1)then
write(9,901),str_quality_1(m8),bondq(m8),str_quality_2(m8),'|',1,1,(astr(m8,m9),m9=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
!write(9,909),quality_fac(m8),1,':',niao,(astr(m8,m9),m9=1,nae)
!endif
tqlty=tqlty+str_quality_1(m8)
bqlty=bqlty+bondq(m8)
sqlty=sqlty+str_quality_2(m8)
enddo
write(9,914)'qualities:','intra_bond','sym_break','nn_bond'
write(9,915),tqlty,sqlty,bqlty

900 format(I3,x,I3,x,I3,x,a,3x,25I4)
901 format(I3,x,I3,x,I3,x,a,3x,I1,a,I2,x,25I4)
909 format(I3,x,I1,a,I1,x,25I4)
914 format(a,x,a,x,a,x,a)
915 format(15x,I3,7x,I3,7x,I3)
goto 374
endif
!**********************************************************************************************************
!!!!!!!!!!!!!!!! If number of permisible structures and available structures are same END !!!!!!!!!!!!!!!
!**********************************************************************************************************

cnt=0
do i5=1,10000
str_cnt(i5)=0
enddo

do i1=1,m4
do m8=1,10000
do m9=1,nae
str1(m8,m9)=0
enddo
enddo

if(repl.ne.0)then
do i4=1,repl
rep=0
do i3=1,nae
if(repl_struc(i4,i3).eq.astr(i1,i3))then
rep=rep+1
endif
enddo 
if(rep.eq.nae)goto 373
enddo
endif


str1_cnt=0

do i5=1,cnt
if(str_cnt(i5).eq.i1)goto 373
enddo

cnt=cnt+1
str1_cnt=str1_cnt+1

do i3=1,nae
str1(str1_cnt,i3)=astr(i1,i3)
enddo 

!print*,'sourav r',str1_cnt,i1
symqq(str1_cnt)=symq(i1)
qulsymm(str1_cnt)=qulsym(i1)
!print*,'symqq(str1_cnt),qulsymm(str1_cnt)',symqq(str1_cnt),qulsymm(str1_cnt)
!print*,symqq(str1_cnt)

lp_cnt=0

do i3=1,nl*2,2
lp_cnt=lp_cnt+1
lps(lp_cnt)=astr(i1,i3)
enddo

!write(*,305),str1_cnt,str1_cnt,(str1(str1_cnt,m9),m9=1,nae)
str_cnt(cnt)=i1

do i2=i1,m4

if(repl.ne.0)then
do i4=1,repl
rep=0
!print*,(repl_struc(i4,i3),i3=1,nae)
do i3=1,nae
if(repl_struc(i4,i3).eq.astr(i2,i3))then
rep=rep+1
endif
enddo 
if(rep.eq.nae)goto 372
enddo
endif

do i5=1,cnt
if(str_cnt(i5).eq.i2)goto 372
enddo
!write(*,305),i2,(astr(i2,m9),m9=1,nae)


if(nl.ne.0)then
ii=0
do i6=1,nl*2,2
do i7=1,nl*2,2
if(astr(i1,i6).eq.astr(i2,i7))then
ii=ii+1
endif
enddo
enddo
if(ii.ne.nl)goto 372
endif


do i3=nl*2+1,nae
do i4=nl*2+1,nae
if(astr(i1,i3).eq.astr(i2,i4))goto 371
enddo
goto 372
371 enddo

cnt=cnt+1
str1_cnt=str1_cnt+1

do i3=1,nae
str1(str1_cnt,i3)=astr(i2,i3)
enddo

symqq(str1_cnt)=symq(i2)
qulsymm(str1_cnt)=qulsym(i2)
!print*,'symqq(str1_cnt),qulsymm(str1_cnt)',symqq(str1_cnt),qulsymm(str1_cnt)

str_cnt(cnt)=i2

elporb=nae-nl*2
total_str=dfactorial(elporb-1-mod(elporb,2))
!print*,'total_str',total_str,elporb,nl,nae
!print*,'cnt',nl,total_str,elporb
if(flg_ion.eq.1)then
if(total_str.eq.cnt)goto 375
endif
if(flg_cov.eq.1)then
if(total_str*elporb.eq.cnt)goto 375
endif
372 enddo 

!elporb=nae-nl*2
375 call wigner(elporb,wig2)
write(7,*),wig2,' structures are permisible among below ',str1_cnt
if(nl.ne.0)write(7,*)'                 lone pair =',(lps(i3),i3=1,lp_cnt)

call rumer_structures(nl,str1,str1_cnt,rumer,rumer_rad)
call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq)
!call vector_rep(nl,str1,str1_cnt,fvec)

!call symmetry_cal(nl,str1,str1_cnt,symsc,symq,nssym)

!n=0
!do ijk=nssym,1,-1
do m8=1,str1_cnt
!if(symq(m8).eq.ijk)then
!n=n+1
!write(7,*),rumer(m8),rumer_rad(m8),rumer(m8)*rumer_rad(m8)
if(rumer(m8)*rumer_rad(m8).eq.1)rumstr='R'
if(rumer(m8)*rumer_rad(m8).eq.0)rumstr='-'
write(7,301),'str',m8,')','[',str_quality_1(m8),',',str_quality_2(m8),']','{',quality_fac(m8),'}',&
rumstr,(str1(m8,m9),m9=1,nae)
!write(7,305)(fvec(m8,m9),m9=1,ndet)
!write(*,305)(str1(m8,m9),m9=1,nae)
!endif
enddo
!write(7,*)
!enddo

!call mat_ind(nl,str1,str1_cnt)

!***********************************************************************
!!!!!!!!!!!!!!!! Rumer Structures selection  Start !!!!!!!!!!!!!!!!!!!!!
!***********************************************************************

if(flg1.eq.1)then

call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq)

!do i=1,str1_cnt
!write(*,231),(str1(i,j),j=1,nae),str_quality_1(i),str_quality_2(i),bondq(i),quality_fac(i)
!enddo


!do i=1,str1_cnt
!write(*,231),(str2(i,j),j=1,nae),rumer(i),rumer_rad(i)
!enddo

call qult_str_arrange(nl,str1,str1_cnt,quality_fac,str2,q_fac)
call rumer_structures(nl,str2,str1_cnt,rumer,rumer_rad)
!call vector_rep(nl,str2,str1_cnt,fvec)
call write_rumer_xmi(nl,str2,str1_cnt,rumer,rumer_rad,quality_fac)
endif
print*,'sourav1'
!***********************************************************************
!!!!!!!!!!!!!!!! Rumer Structures selection Ends !!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!***********************************************************************
!!!!!!!!!!!!!!!! Ionic Structures selection for CHEM. QUAL. Starts !!!!!!!!!!!!!!!!!!
!***********************************************************************

!print*,'******** flg_ion,flg_cov',flg_ion,flg_cov
if(flg1.eq.0.and.flg_ion.eq.1.and.flg_cov.eq.0.and.nfset.ne.0)then
call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq)

!do i=1,str1_cnt
!write(*,231),(str1(i,j),j=1,nae)
!enddo

!call vector_rep(nl,str1,str1_cnt,fvec)
call qult_str_arrange(nl,str1,str1_cnt,quality_fac,str2,q_fac)
!print*,'Ionic Structures selection'
call vector_rep(nl,str2,str1_cnt,fvec)
if(ovopt.eq.1)call overlap(nl,str2,str1_cnt)
call main_bond_cal(nl,str2,str1_cnt,mbondq)
call write_xmi_new_2(nl,wig2,str2,str1_cnt,q_fac,quality_fac,fvec)
!call write_xmi_new(nl,wig2,str2,str1_cnt,q_fac,tqlty)
!call write_xmi2(nl,wig2,str2,str1_cnt,q_fac,tqlty)
!print*,'tqlty',tqlty

endif

!***********************************************************************
!!!!!!!!!!!!!!!! Ionic Structures selection for CHEM. QUAL. Ends !!!!!!!!!!!!!!!!!!
!***********************************************************************

!***********************************************************************
!!!!!!!!!!!!!!!! Covalent Structures selection Starts !!!!!!!!!!!!!!!!!!
!***********************************************************************

!print*,' flg_ion,flg_cov',flg_ion,flg_cov
if(flg1.eq.0.and.flg_cov.eq.1.and.flg_ion.eq.0)then
!print*,'sourav1'

if(symm.eq.1) then
call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq)


!call sym_qual(str1,str1_cnt,nl,qulsym,tnqs)
!print*,'sourav'tnqs

if(sig_sym_flg.eq.1)call symmetry_cal_sig(nl,str1,str1_cnt,symsc,symq,nssym)
if(sig_sym_flg.ne.1)call symmetry_cal_pi(nl,str1,str1_cnt,symsc,symq,nssym)
!do i=1,str1_cnt
!write(*,231),i,bondq(i),quality_fac(i),qulsymm(i),symq(i),(str1(i,j),j=1,nae)
!enddo
!stop
!call sigma_symm(str1,str1_cnt,nl,sigsym,tnqs_sig)
!print*,'sourav'tnqs
!call vector_rep(nl,str1,str1_cnt,fvec)
call qult_str_arrange(nl,str1,str1_cnt,quality_fac,str2,q_fac)
!call symm_str_arrange(str1,str1_cnt,quality_fac,str2,symqq,nssym,qulsymm,tnqs,q_fac)
!print*,'Covalent Structures selection symm'
do i=1,str1_cnt
write(*,231),(str2(i,j),j=1,nae)
!print*,'str2',quality_fac(i),bondq(i),q_fac(i)
enddo

call vector_rep(nl,str2,str1_cnt,fvec)
if(ovopt.eq.1)call overlap(nl,str2,str1_cnt)
!call write_symm_xmi(nl,wig2,str2,str1_cnt,q_fac)
call write_symm_xmi_new(nl,wig2,str2,str1_cnt,q_fac)
!call write_xmi_new_2(nl,wig2,str2,str1_cnt,q_fac,quality_fac,fvec)

endif

if(symm.eq.0)then
if(nfset.ne.4)then
!print*,'sourav3'
call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq)

!call vector_rep(nl,str1,str1_cnt,fvec)
call qult_str_arrange(nl,str1,str1_cnt,quality_fac,str2,q_fac)

!do i=1,str1_cnt
!write(*,231),(str2(i,j),j=1,nae),q_fac(i),str_quality_1(i),str_quality_2(i),bondq(i)
!print*,'str2',quality_fac(i),bondq(i),q_fac(i)
!enddo
!call charge_cntr(str2,nl,str1_cnt)
!print*,'Covalent Structures selection non_symm'
call vector_rep(nl,str2,str1_cnt,fvec)
!print*,'vector_rep_in'
!call vector_rep(nl,str1,str1_cnt,fvec)
!print*,'vector_rep_out'
!call write_xmi_new_2(nl,wig2,str1,str1_cnt,q_fac1,quality_fac,fvec)
print*,'ovopt',ovopt,niao
if(ovopt.eq.1)call overlap(nl,str2,str1_cnt)
call main_bond_cal(nl,str2,str1_cnt,mbondq)
print*,'niao:str_sel',niao
call write_xmi_new_2(nl,wig2,str2,str1_cnt,q_fac,quality_fac,fvec)

endif

if(nfset.eq.4)then
print*,'str1_cnt',str1_cnt
!stop
call quality_factor(nl,str1,str1_cnt,quality_fac,str_quality_1,str_quality_2,bondq)

!call vector_rep(nl,str1,str1_cnt,fvec)
call qult_str_arrange(nl,str1,str1_cnt,quality_fac,str2,q_fac)

!do i=1,str1_cnt
!write(*,231),(str2(i,j),j=1,nae),q_fac(i),str_quality_1(i),str_quality_2(i),bondq(i)
!print*,'str2',quality_fac(i),bondq(i),q_fac(i)
!enddo

!print*,'Covalent Structures selection non_symm'
!call vector_rep(nl,str2,str1_cnt,fvec)
!print*,'vector_rep_in'
!call vector_rep(nl,str1,str1_cnt,fvec)
!print*,'vector_rep_out'
!call write_xmi_new_2(nl,wig2,str1,str1_cnt,q_fac1,quality_fac,fvec)
print*,'ovopt',ovopt,niao
if(ovopt.eq.1)call overlap(nl,str2,str1_cnt)
call main_bond_cal(nl,str2,str1_cnt,mbondq)
call check_str_bond(nl,wig2,str2,str1_cnt)
print*,'niao:str_sel',niao
call eq_dstr_set(str1_cnt,nl,wig2,str2)
endif
endif

endif

!***********************************************************************
!!!!!!!!!!!!!!!! Covalent Structures selection Ends !!!!!!!!!!!!!!!!!!
!***********************************************************************
tqlty1=tqlty1+tqlty
bqlty1=bqlty1+bqlty
sqlty1=sqlty1+sqlty
!print*,'rajat',tqlty1,tqlty
373 enddo
231 format(30I3)

if(nfset.eq.1.or.flg1.eq.1) then
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
write(9,910)'qualities:',' intra_bond','sym_break','nn_bond'
write(9,911),tqlty1,sqlty1,bqlty1
!write(9,*),'Total Quality = ',tqlty1
write(9,*),'    '
endif

910 format(a,x,a,x,a,x,a)
911 format(15x,I3,7x,I5,7x,I3)

print*,'exit str_selection'
374 return
end subroutine str_selection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mat_ind(nl,numstr,totstr,strno,Ifail,det_inv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use commondat
use commondat1
implicit none

!common/str/str5
common/fail/faiil
!common/sdt/strdet,detmnt,det_sign
integer::nl,str1(10000,15),Ifail,totstr,tndet,vec(10000,1000),numstr,strno(1000),detmnt1(10000,15)&
,i,j,k,l,m,n,i1,comdet(10000),mat(10000,10000),sdet(10000),mat_sign(10000,10000),S00(10000,10000)&
,det_sign1(10000),faiil,det_inv
!,strdet(100000),detmnt(100000,15),det_sign(100000)
integer::col(1000)
DOUBLE PRECISION :: SPL0,SPL,D1(100)
character(len=3)::ind

!write(9,*),'enter mat_ind'
!do l=1,totstr*ndet
!print*,'strdet(l)',l,'<',strdet(l)
!enddo
!totstr=6
!call vector_rep(nl,str1,totstr,vec,detmnt,strdet,det_sign,tndet)
tndet=0
!write(9,*),'totstr,numstr',totstr,numstr
!write(9,103),'mat_ind_strno',ndet*totstr,(strno(i),i=1,numstr)

!print*,'numstr',numstr,totstr,ndet
do i=1,numstr
j=strno(i)
!write(9,*),'structures',i,j,numstr
do l=1,totstr*ndet
!write(9,*),'lll',l
if(strdet(l).ne.j)goto 200
!write(9,*),'strdet(l)',strdet(l),j,totstr*ndet
tndet=tndet+1
do m=1,nae
detmnt1(tndet,m)=detmnt(l,m)
enddo
!write(9,103),'detmnt(l,m)',(detmnt(l,m),m=1,nae),det_sign(l)
det_sign1(tndet)=det_sign(l)
200 enddo
enddo
!do i=1,tndet
!write(9,103),'detmnt1***',(detmnt1(i,m),m=1,nae)
!enddo
!print*,'totstr*ndet,tndet',totstr,ndet,tndet

l=0
n=0
do i=1,tndet
!print*,'det',i,(detmnt(i,j),j=1,nae)
!print*,'lll',l
do i1=1,l
if(comdet(i1).eq.i)goto 101
enddo
n=n+1
m=0
do j=i,tndet
do k=1,nae
if(detmnt1(i,k).ne.detmnt1(j,k))goto 100
enddo
l=l+1
m=m+1
mat(n,m)=strdet(j)
mat_sign(n,m)=det_sign1(j)
comdet(l)=j
!print*,'iii',l,comdet(l),n,m,mat(n,m),mat_sign(n,m)
100 enddo
sdet(n)=m
!print*,'sdet(n)',n,sdet(n)
101 enddo

!do i=1,n
!print*,'**',i,sdet(i)
!write(*,102),(mat(i,j),j=1,sdet(i))
!enddo
!do i=1,n
!write(*,102),(mat_sign(i,j),j=1,sdet(i))
!enddo
102 format(50I5)
103 format(a,2x,50I5)

do i=1,totstr
do j=1,totstr
ind_mat(i,j)=0.0
enddo
enddo

do i=1,totstr
ind_mat(i,i)=ndet
!write(9,*),'mat_ind::ndet',ndet
enddo

do i=1,numstr
!print*,i,'****',i
do j=i+1,numstr
!print*,'j',j
do k=1,n
!print*,'k',k
do l=1,sdet(k)
!print*,'l,sdet(k),mat(k,l),i',l,sdet(k),mat(k,l),i
if(mat(k,l).eq.i)then
do m=1,sdet(k)
!print*,'m,sdet(k),mat(k,m),j',m,sdet(k),mat(k,m),j
if(mat(k,m).eq.j)then
ind_mat(i,j)=ind_mat(i,j)+mat_sign(k,m)*mat_sign(k,l)
!print*,'i,j,k,l,m,ind_mat(i,j)',i,j,k,l,m,ind_mat(i,j),mat_sign(k,m),mat_sign(k,l)
endif
enddo
endif
enddo
enddo
enddo
enddo
!      Do k=1,n
!print*,'sdet(k)',sdet(k)
!       Do i=1,sdet(k)
!        IIs=Nsd(Is,Jdet)
!        Spl0=Dble(mat_sign(k,i))
!print*,'Spl0',k,i,Spl0
!        Do j=i,sdet(k)
!         JJs=Nsd(Js,Jdet)
!         Spl=Spl0*Dble(mat_sign(k,j))
!         S00(j,i)=S00(j,i)+Spl
!        Enddo
!       Enddo
!      Enddo
      Do i=1,numstr
       Do j=i+1,numstr
        ind_mat(j,i)=ind_mat(i,j)
       Enddo
      Enddo

!Do i=1,numstr
!write(*,102),(ind_mat(i,j),j=1,numstr)
!
!!write(*,102),(S00(i,j),j=1,totstr)
!Enddo


!call MatLDR('ind',col,numstr,D1)
!det_inv=1.0
!do i=1,numstr
!if (dabs(D1(i)).lt.1.D-10)then
!Ifail=1
!goto 111
!endif
!
!det_inv=det_inv*D1(i)
!enddo
!print*,'mat_ind',(D1(i),i=1,numstr),'det_inv=',det_inv
!
!if (dabs(det_inv).lt.1.D-10)then
!Ifail=1
!else
!Ifail=0
!endif

call Invmat(numstr,Ifail)
111 faiil=Ifail
!if(numstr.eq.14.and.Ifail.eq.0)then
!print*,(strno(i),i=1,numstr)
!if(Ifail.eq.0)then
!write(9,*),'Ifail:invmat',Ifail
!do i=1,numstr
!write(9,103),'mat_ind_str',(str5(strno(i),j),j=1,nae)
!enddo
!Do i=1,numstr
!write(9,102),(ind_mat(i,j),j=1,numstr)
!!!!!
!!!!!!write(*,102),(S00(i,j),j=1,totstr)
!Enddo
!endif
!endif
!if(Ifail.eq.0.and.numstr.eq.5)serial=serial+1
!print*,'serial=',serial
!print*,'exit mat_ind'

!if(numstr.eq.7)stop
return
end subroutine mat_ind

!---------------------------------------------
      SUBROUTINE Invmat(Ndim,Ifail)
!---------------------------------------------
!
!
!
use commondat
use commondat1
      Implicit DOUBLE PRECISION(A-H,O-Z)
      Dimension A(1000,1000),B(1000,1000)




!print*,'enter Invmat'
do i=1,ndim
do j=1,ndim
A(i,j)=ind_mat(i,j)
!print*,A(i,j)
enddo
enddo
!do i=1,ndim
!write(9,111),'ind_mat',(int(A(i,j)),j=1,ndim)
!enddo
111 format(a,30I5)

      N=Ndim
      Ifail = 0
      Xdet = 1.D0
      Do I = 1,N
       Do J = 1,N
        B(J,I) = 0.D0
       Enddo
       B(I,I) = 1.D0
      Enddo
!do i=1,ndim
!print*,'A(i,j)',(A(i,j),j=1,ndim)
!enddo

      Do I = 1,N-1
!--Find the Max
       Amax = Dabs(A(I,I))
!print*,'AmaxII',I,Amax
       Nmax = I

       Do J = I+1,N
        Aij = Dabs(A(I,J))
!print*,'Aij',I,J,Aij
        If(Aij.Gt.Amax) Then
         Amax = Aij
         Nmax = J
        Endif
       Enddo

!print*,'Amaxij',Amax,Nmax
       If(Amax.Lt.1.D-11) Then
!       If(Amax.Lt.1.D-1) Then
        Ifail = 1
!print*,'Amax',Amax
        Xdet = 0.D0
        Return
       Endif

       If(Nmax.Ne.I) Then
        Xdet = -Xdet
        Do J = 1,N
         Atmp = A(J,I)
         A(J,I) = A(J,Nmax)
         A(J,Nmax) = Atmp
         Atmp = B(J,I)
         B(J,I) = B(J,Nmax)
         B(J,Nmax) = Atmp
        Enddo
       Endif

!-- Vanish Lower
       Do J = I+1,N
        Atmp = A(I,J)/A(I,I)
        Do K = 1,N
!print*,'A(K,J)',I,J,K,A(K,J),Atmp,A(K,I)
         A(K,J) = A(K,J)-Atmp*A(K,I)
         B(K,J) = B(K,J)-Atmp*B(K,I)
!print*,'A(K,J)**',I,J,K,A(K,J),Atmp,A(K,I)
!do i1=1,ndim
!print*,'A(i,j)',(A(i1,j1),j1=1,ndim)
!enddo
        Enddo
       Enddo
        Do K = 1,N
!print*,'A(K,J)',k,(A(K,J),J=1,N)
        Enddo
      Enddo
!do i=1,ndim
!print*,'A(i,j)',(A(i,j),j=1,ndim)
!enddo

      If(Dabs(A(N,N)).Lt.1.D-11) Then
!      If(Dabs(A(N,N)).lt.1.D-1) Then
       Ifail = 1
!print*,'Dabs(A(N,N))',N,Dabs(A(N,N))
       Xdet = 0.D0
       Return
      Endif

      Do I = 1,N
       Xdet = Xdet*A(I,I)
      Enddo

!-- Vanish Upper
      Do I = N,2,-1
       Do J = I-1,1,-1
        Atmp = A(I,J)/A(I,I)
        Do K = 1,N
         A(K,J) = A(K,J)-Atmp*A(K,I)
         B(K,J) = B(K,J)-Atmp*B(K,I)
!print*,'A(K,J)****',I,J,K,A(K,J),Atmp,A(K,I)
        Enddo
       Enddo
      Enddo

!do i=1,ndim
!print*,'A(i,j)',(A(i,j),j=1,ndim)
!enddo
!-- Normalized
      Do I = 1,N
       Atmp = A(I,I)
       Do J = 1,N
        B(J,I) = B(J,I)/Atmp
       Enddo
      Enddo
      End

!--------------------------------------------
      SUBROUTINE MatLDR1(name,col1,N,D)
!--------------------------------------------
!
!     A = U * D * VT**T
!     by two-sides Jacobi's method
!     U,VT::: left and right eigen value of A
!     right eigenvectors are stored in transpose form
!     D::: diagonal matrix stored in one dimension array
!     N::: the length of row/column of matrix A
!     Nlt::: is the Nulity of matrix A, which is equal to the number of zero eigenvalues 
!     Up::: upper main tridiagonal
!     Low::: lowwer main tridiagonal
!     Amax::: maximum off-diagonal matrix elements
!
      USE commondat1
      Implicit none
      Double Precision A(N,N),U(N,N),D(N),VT(N,N)
      Double Precision Amax,Anzo,Aii,Ajj,Aij,Aji,Aod,Adif,tmp
      integer N,Nlt,col1(2,1000)
      integer I,J,K,L,Iplus,faiil
      integer iter,Itmax,Dotri
      Double Precision Zer,One,Eps,ovlp
      INTEGER IW,IER
      character(len=3)::name
common/fail/faiil

      Data Zer,One,Eps/0.D0,1.D0,1.0D-20/

!print*,'name_gandu',name
if(name.eq.'orb')then
do i=1,N
do j=1,N
A(i,j)=orb_ovlp_mat1(col1(1,i),col1(2,j))
enddo
enddo
endif

!if(N.eq.5.and.faiil.eq.0)then
!print*,'col(i)_Mat',(col1(1,i),i=1,N),'col(j)',(col1(2,j),j=1,N)
!print*,'faiil',faiil
do i=1,N
print*,'orb_ovlp_mat1_MatL',(orb_ovlp_mat1(col1(1,i),col1(2,j)),j=1,N)
enddo
!endif
!
!do i=1,N
!print*,'A(i,j)',(A(i,j),j=1,N)
!enddo

!
!if(name.eq.'ind')then
!do i=1,N
!do j=1,N
!A(i,j)=ind_mat(i,j)
!enddo
!print*,'MatLDR:A',(A(i,j),j=1,N)
!enddo
!write(*,102),(ind_mat(i,j),j=1,numstr)
!endif
!      IW  = IVBF(1)
!      IER = IVBF(3)
!     write(0,*)"NDIM = ",Ndim
!     write(0,*)"Size = ",Size(VT,1),Size(VT,2),Size(VT)
!do i=1,N
!print*,'MatLDR:A',(A(i,j),j=1,N)
!enddo

!print*,'N',N,Zer

      do i = 1, N
        D(i) = Zer
        do j = 1, N
          U(j,i) = Zer
          VT(j,i) = Zer
        enddo
        U(i,i) = One
        VT(i,i) = One
      enddo
      Nlt = 0
      if(N.le.0) Return
      if(N.eq.1) then
          D(1) = A(1,1)
          if(dabs(D(1)) .lt. 1.d-12 ) Nlt = 1
          Return
      endif
      itmax = 200
      iter = 0
      Amax = One
      do while(Amax .gt. Eps)
!print*,'Eps',Eps,Amax
          iter = iter + 1
          Amax = Zer
          do i = 2, N
          do j = 1, i - 1
            Aii = A(i,i)
            Ajj = A(j,j)
            Aij = A(i,j)
            Aji = A(j,i)
            if(dabs(Aij) .gt. Amax) Amax = dabs(Aij)
            if(dabs(Aji) .gt. Amax) Amax = dabs(Aji)
            if(dabs(Aij) .lt. Eps .and. dabs(Aji) .lt. Eps) cycle
            Call TwoSideJacobi(i,j,Aii,Ajj,Aij,Aji,U,VT,A,N)
!print*,'Aij',i,j,Aii,Aij,Ajj,Aji
          end do
          end do
      end do
!print*,'Aij***',i,j,Aii,Aij,Ajj,Aji
!do i=1,3
!print*,'U',(U(i,j),j=1,3)
!enddo
!do i=1,3
!print*,'VT',(VT(i,j),j=1,3)
!enddo

      do i = 1,N
        D(i) = A(i,i)
      End do
      if((iter / N) .eq. itmax) then
        write(*,*)'fail in MatLDR'
        write(*,*)'Average number of sweep:::',iter / N
        Stop
      endif

!     bubble sort the eigenvalues and eigenvectors
      Do I = 1, N-1
        tmp = dabs(D(I))
        K = I
        Do J = I+1, N
          If(dabs(D(J)) .GT. tmp) then
            tmp = dabs(D(J))
            K = J
          Endif
        Enddo
        If( K .NE. I) Then
          tmp = D(I)
          D(I) = D(K)
          D(K) = tmp
          Do J = 1, N
            tmp = U(J,I)
            U(J,I) = U(J,K)
            U(J,K) = tmp
            tmp = VT(J,I)
            VT(J,I) = VT(J,K)
            VT(J,K) = tmp
          Enddo
        Endif
      Enddo
      Do I=N,1,-1
        If(Dabs(D(I)).LT.1.d-12) Then
          Nlt=Nlt+1
        Else  
          Exit
        Endif
      Enddo
      Call MatTran(VT,N)
ovlp=1.0
do i=1,N
ovlp=ovlp*D(i)
enddo
!print*,'MatDLR:D',(D(i),i=1,N),ovlp
      Return
      End

!--------------------------------------------
      SUBROUTINE MatLDR(name,col,N,D)
!--------------------------------------------
!
!     A = U * D * VT**T
!     by two-sides Jacobi's method
!     U,VT::: left and right eigen value of A
!     right eigenvectors are stored in transpose form
!     D::: diagonal matrix stored in one dimension array
!     N::: the length of row/column of matrix A
!     Nlt::: is the Nulity of matrix A, which is equal to the number of zero eigenvalues 
!     Up::: upper main tridiagonal
!     Low::: lowwer main tridiagonal
!     Amax::: maximum off-diagonal matrix elements
!
      USE commondat1
      Implicit none
      Double Precision A(N,N),U(N,N),D(N),VT(N,N),bond_vec_mat(150,150)
      Double Precision Amax,Anzo,Aii,Ajj,Aij,Aji,Aod,Adif,tmp
      integer N,Nlt,col(1000)
      integer I,J,K,L,Iplus,faiil
      integer iter,Itmax,Dotri
      Double Precision Zer,One,Eps,ovlp
      INTEGER IW,IER
      character(len=3)::name
common/fail/faiil
common /bnd_mat/bond_vec_mat

      Data Zer,One,Eps/0.D0,1.D0,1.0D-20/

!print*,'name_gandu',name
!if(name.eq.'orb')then
!do i=1,N
!do j=1,N
!A(i,j)=orb_ovlp_mat1(col1(1,i),col1(2,j))
!enddo
!enddo
!endif
!
!print*,'col(i)',(col1(1,i),i=1,N),'col(j)',(col1(2,j),j=1,N)

if(name.eq.'str')then
do i=1,N
do j=1,N
A(i,j)=ovlp_mat_norm(col(i),col(j))
enddo
enddo
endif

if(name.eq.'ind')then
do i=1,N
do j=1,N
A(i,j)=ind_mat(i,j)
enddo
!print*,'MatLDR:A',(A(i,j),j=1,N)
enddo
!write(*,102),(ind_mat(i,j),j=1,numstr)
endif
if(name.eq.'bnv')then
do i=1,N
do j=1,N
A(i,j)=bond_vec_mat(i,j)
enddo
enddo
endif
!if(N.eq.14.and.name.eq.'str'.and.faiil.eq.0)then
!print*,'col(i)',(col(i),i=1,N),'col(j)',(col(j),j=1,N)
!do i=1,N
!print*,'orb_ovlp_mat1_MatL',(ovlp_mat_norm(col(i),col(j)),j=1,N)
!enddo
!endif
!      IW  = IVBF(1)
!      IER = IVBF(3)
!     write(0,*)"NDIM = ",Ndim
!     write(0,*)"Size = ",Size(VT,1),Size(VT,2),Size(VT)
!do i=1,N
!print*,'MatLDR:A',(A(i,j),j=1,N)
!enddo

!print*,'N',N,Zer

      do i = 1, N
        D(i) = Zer
        do j = 1, N
          U(j,i) = Zer
          VT(j,i) = Zer
        enddo
        U(i,i) = One
        VT(i,i) = One
      enddo
      Nlt = 0
      if(N.le.0) Return
      if(N.eq.1) then
          D(1) = A(1,1)
          if(dabs(D(1)) .lt. 1.d-12 ) Nlt = 1
          Return
      endif
      itmax = 200
      iter = 0
      Amax = One
      do while(Amax .gt. Eps)
!print*,'Eps',Eps,Amax
          iter = iter + 1
          Amax = Zer
          do i = 2, N
          do j = 1, i - 1
            Aii = A(i,i)
            Ajj = A(j,j)
            Aij = A(i,j)
            Aji = A(j,i)
            if(dabs(Aij) .gt. Amax) Amax = dabs(Aij)
            if(dabs(Aji) .gt. Amax) Amax = dabs(Aji)
            if(dabs(Aij) .lt. Eps .and. dabs(Aji) .lt. Eps) cycle
            Call TwoSideJacobi(i,j,Aii,Ajj,Aij,Aji,U,VT,A,N)
!print*,'Aij',i,j,Aii,Aij,Ajj,Aji
          end do
          end do
      end do
!print*,'Aij***',i,j,Aii,Aij,Ajj,Aji
!do i=1,3
!print*,'U',(U(i,j),j=1,3)
!enddo
!do i=1,3
!print*,'VT',(VT(i,j),j=1,3)
!enddo

      do i = 1,N
        D(i) = A(i,i)
      End do
      if((iter / N) .eq. itmax) then
        write(*,*)'fail in MatLDR'
        write(*,*)'Average number of sweep:::',iter / N
        Stop
      endif

!     bubble sort the eigenvalues and eigenvectors
      Do I = 1, N-1
        tmp = dabs(D(I))
        K = I
        Do J = I+1, N
          If(dabs(D(J)) .GT. tmp) then
            tmp = dabs(D(J))
            K = J
          Endif
        Enddo
        If( K .NE. I) Then
          tmp = D(I)
          D(I) = D(K)
          D(K) = tmp
          Do J = 1, N
            tmp = U(J,I)
            U(J,I) = U(J,K)
            U(J,K) = tmp
            tmp = VT(J,I)
            VT(J,I) = VT(J,K)
            VT(J,K) = tmp
          Enddo
        Endif
      Enddo
      Do I=N,1,-1
        If(Dabs(D(I)).LT.1.d-12) Then
          Nlt=Nlt+1
        Else  
          Exit
        Endif
      Enddo
      Call MatTran(VT,N)
ovlp=1.0
do i=1,N
ovlp=ovlp*D(i)
enddo
!print*,'MatDLR:D',(D(i),i=1,N)!,1.0-ovlp
      Return
      End
!-------------------------------
      SUBROUTINE MatTran(A,N)
!-------------------------------
!
!
!
        implicit none
        integer N,I,J
        Double Precision A,x
        Dimension A(N,N)

        Do I=1,N-1
          Do J=I+1,N
            x=A(J,I)
            A(J,I)=A(I,J)
            A(I,J)=x
          Enddo
        Enddo

        Return
      End



!----------------------------------------------------------------
      SUBROUTINE TwoSideJacobi(i,j,Aii,Ajj,Aij,Aji,U,VT,A,N)
!----------------------------------------------------------------
!
!     This SUBROUTINE calculate two sides Jacobi rotation to following
!     matrix
!     | Aii Aij |       | Aii~  0  |
!     |         |  -->  |          |
!     | Aji Ajj |       |  0  Ajj~ |
!
      Implicit none
      Double Precision A(N,N),U(N,N),D(N),VT(N,N)
      Double Precision t1,c1,s1
      Double Precision t2,c2,s2
      Double Precision tp,cp,sp
      Double Precision tm,cm,sm
      Double Precision Aii,Ajj,Aij,Aji,tt,Aod,Adif,tmp
      integer N
      integer I,J,k
      Double Precision Zer,One,Eps
      Data Zer,One,Eps/0.D0,1.D0,1.0D-20/
!     write(0,*)"Size 2 = ",Size(VT,1),size(vt,2),size(vt)
!print*,'A11,A12,A21,A22',A(1,1),A(1,2),A(2,1),A(2,2)

!     angle alpha + beta
!     write(0,*)'I J = ',I,J
!     write(0,*)'aii = ',aii
!     write(0,*)'ajj = ',ajj
!     write(0,*)'aij = ',aij
!     write(0,*)'aji = ',aji
      Adif = Aii - Ajj
      Aod = Aij + Aji
      if(dabs(Adif) .lt. dabs(Aod)) then
!       tt is cot(2x)
!       if(tt > Zer) then
!         tp = dsqrt(One + tt*tt) - tt
!       else
!         tp =-dsqrt(One + tt*tt) - tt
!       endif
        tt = Adif / Aod
!       write(0,*)'alpha + beta adif < aod, tt = ',tt
        tp = dsqrt(One + tt*tt)
        if(tt .lt. Zer) tp = -tp
        tp = tp - tt
      else if (dabs(Adif).lt.Eps) then
        tp = Zer
      else
!       tt is tan(2x)
        tt = Aod / Adif
!       write(0,*)'alpha + beta adif >= aod, tt = ',tt
        if(dabs(tt) .gt. 1.d-7) then
          tp = (dsqrt(One + tt*tt) - One) / tt
        else
!         numerical stable:::Taylor expansion
          tp = 0.5d0*tt - 0.125d0*tt*tt + 0.0625d0*tt*tt*tt 
        endif
      endif
      cp = One / dsqrt(One + tp*tp)
      sp = tp * cp
!     angle alpha - beta
      Adif = Aii + Ajj
      Aod = Aji - Aij
      if(dabs(Adif) .lt. dabs(Aod)) then
        tt = Adif / Aod
!       write(0,*)'alpha - beta adif < aod, tt = ',tt
        tm = dsqrt(One + tt*tt)
        if(tt .lt. Zer) tm = -tm
        tm = tm - tt
      else if (dabs(Adif).lt.Eps) then
        tm = Zer
      else
        tt = Aod / Adif
!       write(0,*)'alpha - beta adif >= aod, tt = ',tt
        if(dabs(tt) .gt. 1.d-7) then
          tm = (dsqrt(One + tt*tt) - One) / tt
        else
!         numerical stable:::Taylor expansion
          tm = 0.5d0*tt - 0.125d0*tt*tt + 0.0625d0*tt*tt*tt 
        endif
      endif
      cm = One / dsqrt(One + tm*tm)
      sm = tm * cm
      c1 = cp*cm - sp*sm
      s1 = sp*cm + cp*sm
      c2 = cp*cm + sp*sm
      s2 = sp*cm - cp*sm
!     update U and A
      Do k = 1, N
        tmp=c1*U(k,j)-s1*U(k,i)
        U(k,i)=s1*U(k,j)+c1*U(k,i)
        U(k,j)=tmp
        tmp=c1*A(j,k)-s1*A(i,k)
!print*,A(i,k),A(j,k),s1,c1
        A(i,k)=s1*A(j,k)+c1*A(i,k)
        A(j,k)=tmp
!print*,'A_two',A(i,k),A(j,k),s1,c1,k
      enddo
!     update VT and A
      Do k = 1, N
        tmp=c2*VT(k,j)-s2*VT(k,i)
        VT(k,i)=s2*VT(k,j)+c2*VT(k,i)
        VT(k,j)=tmp
        tmp=c2*A(k,j)-s2*A(k,i)
        A(k,i)=s2*A(k,j)+c2*A(k,i)
        A(k,j)=tmp
!print*,'A_two_2',A(k,i),A(k,j)
      enddo

      A(j,i) = Zer
      A(i,j) = Zer

      return
      End
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine write_symm_xmi(nl,strn,str3,ncqs,q_fac2)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!use commondat
!implicit none
!common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
!sigsym,tnqs_sig
!common/str/str5,nstr7
!
!integer::nl,strn,ncqs,tostr,initstr,i,i1,i2,i3,i4,i5,i6,i7,i8,i9,m119,m18,m19,m20,m21,m23,m24,count&
!,qul(100),nqul,j,jj,jjj,fg,flg,ii5,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,x,y,strs&
!,ymd,alstr
!integer::i5up16,i16i7,m16m21,i5up15,i15i7,m15m21,i5up14,i14i7,m14m21,i5up13,i13i7,m13m21,i5up12,i12i7,m12m21,&
!i5up11,i11i7,m11m21,i5up10,i10i7,m10m21,i5up9,i9i7,m9m21,i5up8,i8i7,m8m21,i5up7,i7i7,m7m21,i5up6,i6i7,m6m21,&
!i5up5,i5i7,m5m21,i5up4,i4i7,m4m21,i5up3,i3i7,m3m21,i5up2,i2i7,m2m21,i5up1,i1i7,m1m21,i5up17,i17i7,m17m21
!integer::str3(15000,20),q_fac2(15000),finalvec(15000),strset(1000),col(1000),sigsym(15000),tnqs_sig,&
!ffvec2(15000,1000),bondq(15000),bondq4(15000),nqset(15000),str5(2000,20),nstr7,&
!tndet,totstr,Ifail,indpnt,strno(1000),str_quality_1(15000),str_quality_2(15000),ttqlty0,ttqlty&
!,tqlty,bqlty,sqlty,hqlty,tnqs,nssym,qulsym(15000),symq(15000),set_number,ttqlty1,det_inv,ttqlty2,ttqlty3
!integer::rumer(15000),rumer_rad(15000),quality_fac(15000),rumset,u1,max_set,mns
!!integer,dimension(:),allocatable::qq1,qq2,qq
!integer::qq1(5000),qq2(5000),qq(5000),str2(2000,20),Rid,set_num(100)
!real*8::ovlp
!character(10)::dd,a
!Double Precision::D(1000)
!character(len=100)::outfile
!
!print*,'enter write_symm_xmi'
!mns=0
!u1=1
!max_set=75000
!
!if(nfset.eq.3.or.nfset.eq.5)then
!rumset=0
!call rumer_structures(nl,str3,ncqs,rumer,rumer_rad)
!call write_rumer_xmi(nl,str3,ncqs,rumer,rumer_rad,quality_fac)
!endif
!
!set_number=0
!bqlty=0
!tqlty=0
!sqlty=0
!hqlty=0
!indpnt=2
!!ovlpval=1.0
!if(noq0.gt.strn)then
!ttqlty0=noq0
!else
!ttqlty0=strn+noq0
!endif
!ttqlty1=strn+noq1
!ttqlty2=strn+noq2
!ttqlty3=strn+noq3
!
!do i=1,ncqs
!write(*,231),i,(str3(i,j),j=1,nae),q_fac2(i),str_quality_1(i),str_quality_2(i),bondq(i)
!enddo
!!if(indpnt.eq.1)then
!!call vector_rep(nl,str3,ncqs,ffvec2)
!!endif
!!print*,'sou1',strn,ndet,nl,ncqs
!!do i1=1,ncqs
!!print*,(ffvec2(i1,i),i=1,ndet)
!!enddo
!!call quality_factor(nl,str3,ncqs,q_fac2,q_1,q_2)
!jj=1
!do m19=1,ncqs
!!print*,'q_fac2',q_fac2(m19)
!if(m19.eq.1)qul(1)=q_fac2(1)
!j=jj
!do i=1,j
!if(qul(i).eq.q_fac2(m19))goto 373
!enddo
!jj=jj+1
!qul(i)=q_fac2(m19)
!!print*,qul(i)
!373 enddo
!nqul=jj
!!print*,'ncqs',ncqs
!
!!call date_and_time(dd)
!!read (dd,'(I8)') ymd
!!if(ymd.eq.20230516)stop
!
!do i=1,nqul
!jjj=0
!jj=0
!do m19=1,ncqs
!!print*,qul(i),q_fac2(m19)
!if(qul(i).eq.q_fac2(m19))then
!jjj=jjj+1
!jj=m19
!endif
!enddo
!nqset(i)=jjj
!strset(i)=jj
!print*,nqset(i),strset(i)
!enddo
!
!flg=0
!totstr=0
!i7=0
!m21=0
!i5=ndet
!do i8=1,1000
!finalvec(i8)=0
!enddo
!do m19=1,10000
!do m18=1,15
!str2(m19,m18)=0
!enddo
!enddo
!
!
!i5up1=i5
!i1i7=i7
!m1m21=m21
!
!
!do m1=1,nqul
!
!  do i=i5up1+1,i5
!  finalvec(i)=0
!  enddo
!
!  do m19=i1i7+1,i7
!  strno(m19)=0
!    do i8=1,15
!    str2(m19,i8)=0
!    enddo
!  enddo
!
!i7=i1i7
!totstr=i1i7
!m21=m1m21
!
!jj=0
!
!flg=0
!
!
!if(nqset(m1).gt.strn)goto 701
!
!
! if(m1.eq.1)then
!  do i=1,strset(m1)
!   totstr=totstr+1
!   strno(totstr)=i
!  enddo
! else
!  do i=strset(m1-1)+1,strset(m1)
!   totstr=totstr+1
!   strno(totstr)=i
!  enddo
! endif
!if(totstr.gt.strn)goto 701
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!
!if(Ifail.eq.1)goto 701
! if(totstr.gt.1.and.ovopt.eq.vpt) then
!   call MatLDR('str',strno,totstr,D)
!   ovlp=1.0
!  do i=1,totstr
!   ovlp=ovlp*D(i)
!  enddo
!  if(1.0-ovlp.gt.ovval)goto 701
! endif
!
!
!
! m21=m21+nqset(m1)
! if(m21.gt.strn) goto 701
!
!
! if(m1.eq.1)strs=0
!  if(m1.ne.1)strs=strset(m1-1)
!
! do m19=strs+1,strset(m1)
!   if(q_fac2(m19).ne.qul(m1))goto 601
!   i7=i7+1
! do i3=1,nae
!   str2(i7,i3)=str3(m19,i3)
! enddo
!  col(i7)=m19
!  qq(i7)=q_fac2(m19)
!  qq1(i7)=str_quality_1(m19)
!  qq2(i7)=str_quality_2(m19)
!  bondq4(i7)=bondq(m19)
!  flg=1
!  count=1
! 601 enddo
! if(i7.eq.strn) then
!   ttqlty=0
!   tqlty=0
!   bqlty=0
!   sqlty=0
! do m119=1,i7
!  ttqlty=ttqlty+qq(m119)
! enddo
! if(ttqlty.le.ttqlty0)then
!  mns=mns+1
!  set_number=set_number+1
!  if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!  if(nfset.eq.1.and.set_number.gt.1)then
!   if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!  endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1) goto 379
!goto 701
!endif
!
!!endif
!
!i5up2=i5
!i2i7=i7
!m2m21=m21
!
!do m2=m1+1,nqul
!
!do i=i5up2+1,i5
!finalvec(i)=0
!enddo
!
!do m19=i2i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!
!i7=i2i7
!totstr=i2i7
!m21=m2m21
!
!jj=0
!
!flg=0
!!do m19=strset(m2-1)+1,strset(m2)
!
!if(nqset(m2).gt.strn)goto 702
!
!do i=strset(m2-1)+1,strset(m2)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!if(totstr.gt.strn)goto 702
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!
!if(Ifail.eq.1) goto 702
!!do i=strset(m2-1)+1,strset(m2)
!!strno(totstr)=0
!!totstr=totstr-1
!!enddo
!!goto 702
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 702
!endif
!m21=m21+nqset(m2)
!if(m21.gt.strn) goto 702
!!print*,'sourav_xmi_sym'
!
!
!do m19=strset(m2-1)+1,strset(m2)
!if(q_fac2(m19).ne.qul(m2))goto 602
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=2
!602 enddo
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:','
!!intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 702
!endif
!
!
!i5up3=i5
!i3i7=i7
!m3m21=m21
!
!
!do m3=m2+1,nqul
!
!do i=i5up3+1,i5
!finalvec(i)=0
!enddo
!
!do m19=i3i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!
!i7=i3i7
!totstr=i3i7
!m21=m3m21
!
!jj=0
!
!flg=0
!!do m19=strset(m3-1)+1,strset(m3)
!
!if(nqset(m3).gt.strn)goto 703
!
!do i=strset(m3-1)+1,strset(m3)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!if(totstr.gt.strn)goto 703
!do i=1,totstr
!enddo
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!
!if(Ifail.eq.1) goto 703
!!do i=strset(m3-1)+1,strset(m3)
!!strno(totstr)=0
!!totstr=totstr-1
!!enddo
!!goto 703
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 703
!endif
!
!m21=m21+nqset(m3)
!if(m21.gt.strn) goto 703
!!print*,'sourav_xmi_sym'
!
!!303 enddo
!
!
!do m19=strset(m3-1)+1,strset(m3)
!if(q_fac2(m19).ne.qul(m3))goto 603
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=3
!603 enddo
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(ovopt.eq.vpt) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:','
!!intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 703
!endif
!
!i5up4=i5
!i4i7=i7
!m4m21=m21
!
!do m4=m3+1,nqul
!!print*,'m4m4m4',m4
!do i=i5up4+1,i5
!finalvec(i)=0
!enddo
!do m19=i4i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!i7=i4i7
!totstr=i4i7
!m21=m4m21
!
!jj=0
!flg=0
!!do m19=strset(m4-1)+1,strset(m4)
!
!if(nqset(m4).gt.strn)goto 704
!
!do i=strset(m4-1)+1,strset(m4)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!
!if(totstr.gt.strn)goto 704
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!!print*,'ifail4',ifail
!if(Ifail.eq.1) goto 704
!!strno(totstr)=0
!!totstr=totstr-1
!!goto 304
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 704
!endif
!
!
!m21=m21+nqset(m4)
!if(m21.gt.strn) goto 704
!!print*,'sourav_xmi_sym'
!
!!304 enddo
!
!
!do m19=strset(m4-1)+1,strset(m4)
!if(q_fac2(m19).ne.qul(m4))goto 604
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=4
!604 enddo
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(ovopt.eq.vpt) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:','
!!intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 704
!endif
!
!i5up5=i5
!i5i7=i7
!m5m21=m21
!
!!!!!! loop 5 >>>
!do m5=m4+1,nqul
!!print*,'m5m5m5',m5
!do i=i5up5+1,i5
!finalvec(i)=0
!enddo
!do m19=i5i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!i7=i5i7
!totstr=i5i7
!!print*,'totstr5',totstr
!m21=m5m21
!jj=0
!
!flg=0
!!do m19=strset(m5-1)+1,strset(m5)
!
!if(nqset(m5).gt.strn)goto 705
!
!do i=strset(m5-1)+1,strset(m5)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!
!if(totstr.gt.strn)goto 705
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!if(Ifail.eq.1) goto 705
!
!!strno(totstr)=0
!!totstr=totstr-1
!!goto 305
!!endif
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 705
!endif
!
!
!m21=m21+nqset(m5)
!if(m21.gt.strn) goto 705
!!print*,'sourav_xmi_sym'
!
!!305 enddo
!
!!print*,'jj',jj,m5,nqul
!!print*,'******jj',jj
!
!do m19=strset(m5-1)+1,strset(m5)
!if(q_fac2(m19).ne.qul(m5))goto 605
!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=5
!605 enddo
!!print*,'i7 5',i7,m5
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:','
!!intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 705
!endif
!
!i5up6=i5
!i6i7=i7
!m6m21=m21
!
!do m6=m5+1,nqul
!
!!print*,'m6m6m6',m6
!do i=i5up6+1,i5
!finalvec(i)=0
!enddo
!do m19=i6i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!i7=i6i7
!totstr=i6i7
!!print*,'totstr6',totstr
!m21=m6m21
!
!jj=0
!
!flg=0
!!do m19=strset(m6-1)+1,strset(m6)
!
!if(nqset(m6).gt.strn)goto 706
!
!do i=strset(m6-1)+1,strset(m6)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!
!if(totstr.gt.strn)goto 706
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!!print*,'ifail6',ifail
!if(Ifail.eq.1)goto 706
!!!do i8=1,nae
!!!str4(totstr,i8)=0
!!!enddo
!!strno(totstr)=0
!!totstr=totstr-1
!!goto 306
!!endif
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 706
!endif
!
!
!m21=m21+nqset(m6)
!if(m21.gt.strn) goto 706
!!print*,'sourav_xmi_sym'
!
!!306 enddo
!
!!print*,'jj',jj,m6,nqul
!do m19=strset(m6-1)+1,strset(m6)
!if(q_fac2(m19).ne.qul(m6))goto 606
!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=6
!606 enddo
!!print*,'i7  6',i7,m6
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0,tqlty,ttqlty1,bqlty,ttqlty2,sqlty,ttqlty3
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'qualities:','
!!intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 706
!endif
!
!i5up7=i5
!i7i7=i7
!m7m21=m21
!
!do m7=m6+1,nqul
!!print*,'m7m7m7',m7
!do i=i5up7+1,i5
!finalvec(i)=0
!enddo
!do m19=i7i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!i7=i7i7
!totstr=i7i7
!!print*,'totstr7',totstr
!m21=m7m21
!
!jj=0
!
!flg=0
!!do m19=strset(m7-1)+1,strset(m7)
!
!if(nqset(m7).gt.strn)goto 707
!do i=strset(m7-1)+1,strset(m7)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!if(totstr.gt.strn)goto 707
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!if(Ifail.eq.1)goto 707
!!!do i8=1,nae
!!!str4(totstr,i8)=0
!!!enddo
!!strno(totstr)=0
!!totstr=totstr-1
!!goto 307
!!endif
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 707
!endif
!
!
!m21=m21+nqset(m7)
!if(m21.gt.strn) goto 707
!!print*,'sourav_xmi_sym'
!
!!307 enddo
!
!!print*,'jj',jj,m7,nqul
!!407 if(jj.eq.nqset(m7))then
!!print*,'******jj',jj
!
!do m19=strset(m7-1)+1,strset(m7)
!if(q_fac2(m19).ne.qul(m7))goto 607
!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=7
!607 enddo
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:','
!!intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 707
!endif
!
!i5up8=i5
!i8i7=i7
!m8m21=m21
!
!do m8=m7+1,nqul
!!print*,'m8m8m8',m8
!do i=i5up8+1,i5
!finalvec(i)=0
!enddo
!do m19=i8i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!i7=i8i7
!totstr=i8i7
!!print*,'totstr8',totstr
!m21=m8m21
!
!jj=0
!
!flg=0
!!do m19=strset(m8-1)+1,strset(m8)
!
!if(nqset(m8).gt.strn)goto 708
!
!do i=strset(m8-1)+1,strset(m8)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!
!if(totstr.gt.strn)goto 708
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!if(Ifail.eq.1)goto 708
!
!!!do i8=1,nae
!!!str4(totstr,i8)=0
!!!enddo
!!strno(totstr)=0
!!totstr=totstr-1
!!goto 308
!!endif
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 708
!endif
!
!
!m21=m21+nqset(m8)
!if(m21.gt.strn) goto 708
!!print*,'sourav_xmi_sym'
!
!!308 enddo
!
!!print*,'jj',jj,m7,nqul
!!408 if(jj.eq.nqset(m8))then
!!print*,'******jj',jj
!
!do m19=strset(m8-1)+1,strset(m8)
!if(q_fac2(m19).ne.qul(m8))goto 608
!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=8
!608 enddo
!!print*,'i7 8',i7,m8
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 708
!endif
!
!i5up9=i5
!i9i7=i7
!m9m21=m21
!
!do m9=m8+1,nqul
!!print*,'m9m9m9',m9
!do i=i5up9+1,i5
!finalvec(i)=0
!enddo
!do m19=i9i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!i7=i9i7
!totstr=i9i7
!!print*,'totstr9',totstr
!m21=m9m21
!
!jj=0
!
!flg=0
!!do m19=strset(m9-1)+1,strset(m9)
!
!if(nqset(m9).gt.strn)goto 709
!
!do i=strset(m9-1)+1,strset(m9)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!
!if(totstr.gt.strn)goto 709
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!!print*,'ifail9',ifail
!if(Ifail.eq.1)goto 709
!!!do i8=1,nae
!!!str4(totstr,i8)=0
!!!enddo
!!strno(totstr)=0
!!totstr=totstr-1
!!goto 309
!!endif
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 709
!endif
!
!
!m21=m21+nqset(m9)
!if(m21.gt.strn) goto 709
!!print*,'sourav_xmi_sym'
!
!!309 enddo
!
!!print*,'jj',jj,m9,nqul
!!409 if(jj.eq.nqset(m9))then
!!print*,'******jj',jj
!
!do m19=strset(m9-1)+1,strset(m9)
!if(q_fac2(m19).ne.qul(m9))goto 609
!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=9
!609 enddo
!!print*,'i7 9',i7,m9
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 709
!endif
!
!i5up10=i5
!i10i7=i7
!m10m21=m21
!
!do m10=m9+1,nqul
!!print*,'m10m10m10',m10
!do i=i5up10+1,i5
!finalvec(i)=0
!enddo
!do m19=i10i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!i7=i10i7
!totstr=i10i7
!!print*,'totstr10',totstr
!m21=m10m21
!
!jj=0
!
!flg=0
!!do m19=strset(m10-1)+1,strset(m10)
!
!if(nqset(m10).gt.strn)goto 710
!
!do i=strset(m10-1)+1,strset(m10)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!
!if(totstr.gt.strn)goto 710
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!if(Ifail.eq.1)goto 710
!!!do i8=1,nae
!!!str4(totstr,i8)=0
!!!enddo
!!strno(totstr)=0
!!totstr=totstr-1
!!goto 310
!!endif
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 710
!endif
!
!
!m21=m21+nqset(m10)
!if(m21.gt.strn) goto 710
!!print*,'sourav_xmi_sym'
!
!!310 enddo
!
!!print*,'jj',jj,m10,nqul
!!410 if(jj.eq.nqset(m10))then
!!print*,'******jj',jj
!
!do m19=strset(m10-1)+1,strset(m10)
!if(q_fac2(m19).ne.qul(m10))goto 610
!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=10
!610 enddo
!!print*,'i7 10',i7,m10
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 710
!endif
!
!i5up11=i5
!i11i7=i7
!m11m21=m21
!
!do m11=m10+1,nqul
!!print*,'m11m11m11',m11
!do i=i5up11+1,i5
!finalvec(i)=0
!enddo
!do m19=i11i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!i7=i11i7
!totstr=i11i7
!!print*,'totstr11',totstr
!m21=m11m21
!
!jj=0
!
!flg=0
!!do m19=strset(m11-1)+1,strset(m11)
!
!if(nqset(m11).gt.strn)goto 711
!
!do i=strset(m11-1)+1,strset(m11)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!
!if(totstr.gt.strn)goto 711
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!if(Ifail.eq.1)goto 711
!!!do i8=1,nae
!!!str4(totstr,i8)=0
!!!enddo
!!strno(totstr)=0
!!totstr=totstr-1
!!goto 311
!!endif
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 711
!endif
!
!
!m21=m21+nqset(m11)
!if(m21.gt.strn) goto 711
!!print*,'sourav_xmi_sym'
!
!!311 enddo
!
!!print*,'jj',jj,m11,nqul
!!411 if(jj.eq.nqset(m11))then
!!print*,'******jj',jj
!
!do m19=strset(m11-1)+1,strset(m11)
!if(q_fac2(m19).ne.qul(m11))goto 611
!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=11
!611 enddo
!!print*,'i7 11',i7,m11
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 711
!endif
!
!i5up12=i5
!i12i7=i7
!m12m21=m21
!
!do m12=m11+1,nqul
!!print*,'m12m12m12',m12
!do i=i5up12+1,i5
!finalvec(i)=0
!enddo
!do m19=i12i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!i7=i12i7
!totstr=i12i7
!!print*,'totstr12',totstr
!m21=m12m21
!
!jj=0
!
!flg=0
!!do m19=strset(m12-1)+1,strset(m12)
!
!if(nqset(m12).gt.strn)goto 712
!
!do i=strset(m12-1)+1,strset(m12)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!
!if(totstr.gt.strn)goto 712
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!if(Ifail.eq.1)goto 712
!!!do i8=1,nae
!!!str4(totstr,i8)=0
!!!enddo
!!strno(totstr)=0
!!totstr=totstr-1
!!goto 312
!!endif
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 712
!endif
!
!
!m21=m21+nqset(m12)
!if(m21.gt.strn) goto 712
!!print*,'sourav_xmi_sym'
!
!
!!312 enddo
!
!!print*,'jj',jj,m12,nqul
!!412 if(jj.eq.nqset(m12))then
!!print*,'******jj',jj
!
!do m19=strset(m12-1)+1,strset(m12)
!if(q_fac2(m19).ne.qul(m12))goto 612
!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=12
!612 enddo
!!print*,'i7 12',i7,m12
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 712
!endif
!
!i5up13=i5
!i13i7=i7
!m13m21=m21
!
!do m13=m12+1,nqul
!
!!print*,'m13m13m13',m13
!do i=i5up13+1,i5
!finalvec(i)=0
!enddo
!do m19=i13i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!i7=i13i7
!totstr=i13i7
!!print*,'totstr13',totstr
!m21=m13m21
!
!jj=0
!
!flg=0
!!do m19=strset(m13-1)+1,strset(m13)
!
!if(nqset(m13).gt.strn)goto 713
!
!do i=strset(m13-1)+1,strset(m13)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!
!if(totstr.gt.strn)goto 713
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!!print*,'ifail13',ifail
!if(Ifail.eq.1)goto 713
!!!do i8=1,nae
!!!str4(totstr,i8)=0
!!!enddo
!!strno(totstr)=0
!!totstr=totstr-1
!!goto 313
!!endif
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 713
!endif
!
!m21=m21+nqset(m13)
!if(m21.gt.strn) goto 713
!!print*,'sourav_xmi_sym'
!
!!313 enddo
!
!!print*,'jj',jj,m13,nqul
!!413 if(jj.eq.nqset(m13))then
!!print*,'******jj',jj
!
!do m19=strset(m13-1)+1,strset(m13)
!if(q_fac2(m19).ne.qul(m13))goto 613
!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=13
!613 enddo
!!print*,'i7 13',i7,m13
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 713
!endif
!
!i5up14=i5
!i14i7=i7
!m14m21=m21
!
!do m14=m13+1,nqul
!!print*,'m14m14m14',m14
!do i=i5up14+1,i5
!finalvec(i)=0
!enddo
!do m19=i14i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!i7=i14i7
!totstr=i14i7
!!print*,'totstr14',totstr
!m21=m14m21
!
!jj=0
!
!flg=0
!!do m19=strset(m14-1)+1,strset(m14)
!
!if(nqset(m14).gt.strn)goto 714
!
!do i=strset(m14-1)+1,strset(m14)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!
!if(totstr.gt.strn)goto 714
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!if(Ifail.eq.1)goto 714
!!!do i8=1,nae
!!!str4(totstr,i8)=0
!!!enddo
!!strno(totstr)=0
!!totstr=totstr-1
!!goto 314
!!endif
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 714
!endif
!
!
!m21=m21+nqset(m14)
!if(m21.gt.strn) goto 714
!!print*,'sourav_xmi_sym'
!
!!314 enddo
!
!!print*,'jj',jj,m14,nqul
!!414 if(jj.eq.nqset(m14))then
!!print*,'******jj',jj
!
!do m19=strset(m14-1)+1,strset(m14)
!if(q_fac2(m19).ne.qul(m14))goto 614
!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=14
!614 enddo
!!print*,'i7 14',i7,m14
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 714
!endif
!
!i5up15=i5
!i15i7=i7
!m15m21=m21
!
!do m15=m14+1,nqul
!!print*,'m15m15m15',m15
!do i=i5up15+1,i5
!finalvec(i)=0
!enddo
!do m19=i15i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!i7=i15i7
!totstr=i15i7
!!print*,'totstr15',totstr
!m21=m15m21
!
!jj=0
!
!flg=0
!!do m19=strset(m15-1)+1,strset(m15)
!
!if(nqset(m15).gt.strn)goto 715
!
!do i=strset(m15-1)+1,strset(m15)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!
!if(totstr.gt.strn)goto 715
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!if(Ifail.eq.1)goto 715
!!!do i8=1,nae
!!!str4(totstr,i8)=0
!!!enddo
!!strno(totstr)=0
!!totstr=totstr-1
!!goto 315
!!endif
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 715
!endif
!
!
!m21=m21+nqset(m15)
!if(m21.gt.strn) goto 715
!!print*,'sourav_xmi_sym'
!
!!315 enddo
!
!!print*,'jj',jj,m15,nqul
!!415 if(jj.eq.nqset(m15))then
!!print*,'******jj',jj
!
!do m19=strset(m15-1)+1,strset(m15)
!if(q_fac2(m19).ne.qul(m15))goto 615
!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=15
!615 enddo
!!print*,'i7 15',i7,m15
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 715
!endif
!
!i5up16=i5
!i16i7=i7
!m16m21=m21
!
!do m16=m15+1,nqul
!!print*,'m16m16m16',m16
!do i=i5up16+1,i5
!finalvec(i)=0
!enddo
!do m19=i16i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!i7=i16i7
!totstr=i16i7
!!print*,'totstr16',totstr
!m21=m16m21
!
!jj=0
!
!flg=0
!!do m19=strset(m16-1)+1,strset(m16)
!
!if(nqset(m16).gt.strn)goto 716
!
!do i=strset(m16-1)+1,strset(m16)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!if(totstr.gt.strn)goto 716
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!if(Ifail.eq.1)goto 716
!!!do i8=1,nae
!!!str4(totstr,i8)=0
!!!enddo
!!strno(totstr)=0
!!totstr=totstr-1
!!goto 316
!!endif
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 716
!endif
!
!
!m21=m21+nqset(m16)
!if(m21.gt.strn) goto 716
!!print*,'sourav_xmi_sym'
!
!!316 enddo
!
!!print*,'jj',jj,m16,nqul
!!416 if(jj.eq.nqset(m16))then
!!print*,'******jj',jj
!
!do m19=strset(m16-1)+1,strset(m16)
!if(q_fac2(m19).ne.qul(m16))goto 616
!!write(9,231),m19,(finalvec(i1),i1=1,i5)
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!count=16
!616 enddo
!!print*,'i7 16',i7,m16
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!goto 716
!endif
!
!
!i5up17=i5
!i17i7=i7
!m17m21=m21
!
!do m17=m16+1,nqul
!!print*,'m17m17m17',m17
!do i=i5up17+1,i5
!finalvec(i)=0
!enddo
!do m19=i17i7+1,i7
!strno(m19)=0
!do i8=1,15
!str2(m19,i8)=0
!enddo
!enddo
!i7=i17i7
!totstr=i17i7
!!print*,'totstr17',totstr
!m21=m17m21
!
!jj=0
!
!flg=0
!!do m19=strset(m17-1)+1,strset(m17)
!
!if(nqset(m17).gt.strn)goto 717
!
!do i=strset(m17-1)+1,strset(m17)
!totstr=totstr+1
!strno(totstr)=i
!enddo
!
!if(totstr.gt.strn)goto 717
!call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!!write(9,*),'ifail',ifail,det_inv
!if(Ifail.eq.1)goto 717
!!strno(totstr)=0
!!totstr=totstr-1
!!goto 317
!!endif
!!endif
!if(totstr.gt.1.and.ovopt.eq.vpt) then
!!print*,'MatLDR_write_xmi',totstr
!call MatLDR('str',strno,totstr,D)
!
!!print*,'DDD',(D(i),i=1,totstr)
!ovlp=1.0
!do i=1,totstr
!ovlp=ovlp*D(i)
!enddo
!!print*,'ovval',1.0-ovlp,ovlp
!if(1.0-ovlp.gt.ovval)goto 717
!endif
!
!
!m21=m21+nqset(m17)
!if(m21.gt.strn) goto 717
!!print*,'sourav_xmi_sym'
!
!!317 enddo
!
!!417 if(jj.eq.nqset(m17))then
!
!do m19=strset(m17-1)+1,strset(m17)
!if(q_fac2(m19).ne.qul(m17))goto 617
!i7=i7+1
!do i3=1,nae
!str2(i7,i3)=str3(m19,i3)
!enddo
!col(i7)=m19
!qq(i7)=q_fac2(m19)
!qq1(i7)=str_quality_1(m19)
!qq2(i7)=str_quality_2(m19)
!bondq4(i7)=bondq(m19)
!flg=1
!617 enddo
!!print*,'i7 16',i7,m17
!if(i7.eq.strn) then
!!if(nfset.gt.0)then
!ttqlty=0
!tqlty=0
!bqlty=0
!sqlty=0
!do m119=1,i7
!ttqlty=ttqlty+qq(m119)
!!tqlty=tqlty+qq1(m119)
!!bqlty=bqlty+bondq4(m119)
!!sqlty=sqlty+qq2(m119)
!enddo
!!endif
!!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0)then
!mns=mns+1
!set_number=set_number+1
!if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
!if(nfset.eq.1.and.set_number.gt.1)then
!if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
!endif
!
!do m19=1,i7
!if(niao.eq.0)then
!write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!enddo
!!if(nfset.eq.0.or.nfset.gt.1) then
!if(nfset.eq.3)then
!Rumwrite=1
!call Rumer_set_id(str2,i7,nl,Rid,set_num)
!endif
!if(ovopt.eq.1) then
!call MatLDR('str',col,i7,D)
!
!ovlp=1.0
!do i=1,i7
!ovlp=ovlp*D(i)
!enddo
!write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
!endif
!!write(9,910)'quality value',ttqlty
!if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
!if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!bqlty=0
!write(9+u1,*),'Set_number=',set_number
!write(9+u1,*),'    '
!endif
!!endif
!!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!!'ovval',1.0-ovlp
!if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
!if(mns.eq.max_set)then
!if(u1.eq.mset-1) goto 379
!close(9+u1)
!mns=0
!u1=u1+1
!if(u1.le.9)then
!write(a,'(I1)')u1
!endif
!if(u1.gt.9.and.u1.lt.100)then
!write(a,'(I2)')u1
!endif
!if(u1.gt.99.and.u1.lt.1000)then
!write(a,'(I3)')u1
!endif
!if(u1.eq.1000)then
!write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
!sets. Please increase the limit'
!goto 379
!endif
!outfile='structure_set_'//trim(a)//trim('.dat')
!open(unit=9+u1,file=outfile,status='unknown')
!endif
!!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!!write(9,911),tqlty,sqlty,bqlty
!!bqlty=0
!!!set_number=set_number+1
!!write(9,*),'Set_number=',set_number
!!write(9,*),'    '
!!endif
!!endif
!!if(nfset.eq.1)goto 379
!endif
!
!717 enddo
!!print*,'count',count
!if(indpnt.eq.2)then
!if (count.eq.1)goto 702
!if (count.eq.2)goto 703
!if (count.eq.3)goto 704
!if (count.eq.4)goto 705
!if (count.eq.5)goto 706
!if (count.eq.6)goto 707
!if (count.eq.7)goto 708
!if (count.eq.8)goto 709
!if (count.eq.9)goto 710
!if (count.eq.10)goto 711
!if (count.eq.11)goto 712
!if (count.eq.12)goto 713
!if (count.eq.13)goto 714
!if (count.eq.14)goto 715
!if (count.eq.15)goto 716
!endif
!!if (count.eq.16)goto 717
!
!if(indpnt.eq.1)then
!if (count.eq.1)goto 701
!if (count.eq.2)goto 702
!if (count.eq.3)goto 703
!if (count.eq.4)goto 704
!if (count.eq.5)goto 705
!if (count.eq.6)goto 706
!if (count.eq.7)goto 707
!if (count.eq.8)goto 708
!if (count.eq.9)goto 709
!if (count.eq.10)goto 710
!if (count.eq.11)goto 711
!if (count.eq.12)goto 712
!if (count.eq.13)goto 713
!if (count.eq.14)goto 714
!if (count.eq.15)goto 715
!if (count.eq.16)goto 716
!endif
!
!716 enddo
!
!715 enddo
!
!714 enddo
!
!713 enddo
!
!712 enddo
!
!711 enddo
!
!710 enddo
!
!709 enddo
!
!708 enddo
!
!707 enddo
!
!706 enddo
!
!705 enddo
!
!704 enddo
!
!703 enddo
!
!702 enddo
!
!701 enddo
!!enddo
!
!
!!379 do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq(m19),(str2(m19,m20),m20=1,nae)
!!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!!write(9,901),qq(m19),bondq4(m19),1,':',niao,(str2(m19,m20),m20=1,nae)
!!!write(9,901),m19,1,':',niao,(str2(m19,m20),m20=1,nae)
!!!write(*,900),m19,niao,(str2(m19,m20),m20=1,nae)
!!!write(9,*)qq(m19)
!!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!!write(9,909),qq(m19),1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!
!
!!do m19=1,i7
!!if(niao.eq.0)then
!!write(9,900),qq1(m19),qq2(m19),bondq4(m19),'|',(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.ne.0)then
!!if(niao.gt.1)then
!!write(9,901),qq1(m19),qq2(m19),bondq4(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!!endif
!!!if(niao.ne.0.and.nnnatom.eq.0)then
!!if(niao.eq.1)then
!!write(9,909),qq1(m19),qq2(m19),bondq4(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!!endif
!!tqlty=tqlty+qq1(m19)
!!bqlty=bqlty+bondq4(m19)
!!sqlty=sqlty+qq2(m19)
!!enddo
!
!!900 format(25I4)
!!901 format(I3,x,I3,x,I1,a,I1,x,25I4)
!!909 format(I3,x,I1,a,I1,x,25I4)
!
!900 format(I3,x,I3,x,I3,x,I3,x,a,x,25I4)
!901 format(I3,x,I3,x,I3,x,I3,x,a,x,I1,a,I3,x,25I4)
!909 format(I3,x,I3,x,I3,x,I3,x,a,x,I3,I3,x,25I4)
!!910 format(a,x,a,x,a,x,a)
!!911 format(15x,I3,7x,I3,7x,I3)
!910 format(a,I3)
!920 format(a,2x,I5,2x,a,2x,100I5)
!911 format(15x,I3,7x,I3,7x,I3)
!912 format(a,3x,F10.3)
!
!close(21)
!
!open(unit=121,file='script10',status='unknown')
!write(121,111)'rm -rf','Rumer_Sets.dat'
!close(121)
!CALL SYSTEM ("chmod +x script10 ")
!CALL SYSTEM ("./script10 ")
!CALL SYSTEM ("rm script10 ")
!111 format(a,x,a)
!
!231 format(50I3)
!!deallocate(str2)
!!deallocate(qq)
!!deallocate(qq1)
!!deallocate(qq2)
!
!379 return
!end subroutine write_symm_xmi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_symm_xmi_new(nl,strn,str3,ncqs,q_fac2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig
common/infosymm/set_number,hqlty,ttqlty0,ttqlty1,ttqlty2,ttqlty3,ncqss,qul,Rid&
,mns,u1,max_set,rumset,nqset,strset,strnn,totstr,strno,incmplt,mincmplt,mincmplt_set,group_num
common/str/str5,nstr7

integer::nl,strn,strnn,ncqs,ncqss,tostr,initstr,i,i1,i2,i3,i4,i5,i6,i7,i8,i9,m119,m18,m19,m20,m21,m23,m24,count&
,qul(100),nqul,j,jj,jjj,fg,flg,ii5,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,x,y,strs
integer::i5up16,i16i7,m16m21,i5up15,i15i7,m15m21,i5up14,i14i7,m14m21,i5up13,i13i7,m13m21,i5up12,i12i7,m12m21,&
i5up11,i11i7,m11m21,i5up10,i10i7,m10m21,i5up9,i9i7,m9m21,i5up8,i8i7,m8m21,i5up7,i7i7,m7m21,i5up6,i6i7,m6m21,&
i5up5,i5i7,m5m21,i5up4,i4i7,m4m21,i5up3,i3i7,m3m21,i5up2,i2i7,m2m21,i5up1,i1i7,m1m21,i5up17,i17i7,m17m21
integer::str3(15000,20),q_fac2(15000),finalvec(15000),strset(1000),col(1000),sigsym(15000),tnqs_sig,&
ffvec2(15000,1000),bondq(15000),bondq4(15000),nqset(15000),str5(2000,20),nstr7,qfac3(15000),group_num(15000)&
,tndet,totstr,Ifail,indpnt,strno(1000),str_quality_1(15000),str_quality_2(15000),ttqlty0,ttqlty&
,tqlty,bqlty,sqlty,hqlty,tnqs,nssym,qulsym(15000),symq(15000),set_number,ttqlty1,det_inv,ttqlty2,ttqlty3
integer::rumer(15000),rumer_rad(15000),quality_fac(15000),rumset,u1,max_set,mns,mincmplt,mincmplt_set(500)
!integer,dimension(:),allocatable::qq1,qq2,qq
integer::qq1(5000),qq2(5000),qq(5000),str2(2000,20),Rid,set_num(100),sf1,sf2,incmplt
real*8::ovlp
Double Precision::D(1000)
character(10)::dd,a
character(len=100)::outfile


print*,'enter_symm_xmi_new'
Rid=0
mincmplt=0
do i=1,500
mincmplt_set(i)=0
enddo
ncqss=ncqs
strnn=strn
incmplt=1
mns=0
u1=1
max_set=75000

if(nfset.eq.3.or.nfset.eq.5)then
rumset=0
call rumer_structures(nl,str3,ncqss,rumer,rumer_rad)
call write_rumer_xmi(nl,str3,ncqss,rumer,rumer_rad,quality_fac)
endif

set_number=0
bqlty=0
tqlty=0
sqlty=0
hqlty=0
indpnt=2
!ovlpval=1.0
if(noq0.gt.strnn)then
ttqlty0=noq0
else
ttqlty0=strnn+noq0
endif
ttqlty1=strnn+noq1
ttqlty2=strnn+noq2
ttqlty3=strnn+noq3

write(*,*)'sl  structures           group_numbers'
do i=1,ncqss
write(*,231)i,(str3(i,j),j=1,nae),q_fac2(i)
!write(*,231),i,(str3(i,j),j=1,nae),q_fac2(i),str_quality_1(i),str_quality_2(i),bondq(i)
group_num(i)=q_fac2(i)
enddo
231 format(30I3)

jj=1
do m19=1,ncqss
!print*,'q_fac2',q_fac2(m19)
if(m19.eq.1)qul(1)=q_fac2(1)
j=jj
do i=1,j
if(qul(i).eq.q_fac2(m19))goto 373
enddo
jj=jj+1
qul(i)=q_fac2(m19)
!print*,qul(i)
373 enddo
nqul=jj
!print*,'nqul',nqul

do i=1,nqul
jjj=0
jj=0
do m19=1,ncqss
!print*,qul(i),q_fac2(m19)
if(qul(i).eq.q_fac2(m19))then
jjj=jjj+1
jj=m19
endif
enddo
nqset(i)=jjj
strset(i)=jj
!print*,'i,nqset(i),strset(i)',i,nqset(i),strset(i)
enddo

flg=0
totstr=0
i7=0
m21=0
i1i7=0
m1m21=0
!i5=ndet
!do i8=1,1000
!finalvec(i8)=0
!enddo
!do m19=1,10000
!do m18=1,15
!str2(m19,m18)=0
!enddo
!enddo


!i5up1=i5
i1i7=i7
m1m21=m21


do m1=1,nqul
!i5=i5up1
!i7=i1i7
!m21=m1m21
i7=0
m21=0

!print*,'loop1'
call write_symm_xmi_1(i7,m21,m1,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 701
if (sf2.eq.1)goto 100

!print*,'i7,m21**',i7,m21
!i5up2=i5
i2i7=i7
m2m21=m21

do m2=m1+1,nqul
!i5=i5up2
i7=i2i7
m21=m2m21

!print*,'loop2'
call write_symm_xmi_1(i7,m21,m2,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 702
if (sf2.eq.1)goto 100

!print*,'22 i7,m21**',i7,m21
!i5up3=i5
i3i7=i7
m3m21=m21

do m3=m2+1,nqul
!i5=i5up3
i7=i3i7
m21=m3m21

!print*,'loop3'
call write_symm_xmi_1(i7,m21,m3,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 703
if (sf2.eq.1)goto 100

!print*,'33 i7,m21**',i7,m21
!stop
!i5up4=i5
i4i7=i7
m4m21=m21

do m4=m3+1,nqul
!i5=i5up4
i7=i4i7
m21=m4m21

!print*,'loop4'
call write_symm_xmi_1(i7,m21,m4,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 704
if (sf2.eq.1)goto 100

!i5up5=i5
i5i7=i7
m5m21=m21

do m5=m4+1,nqul
!i5=i5up5
i7=i5i7
m21=m5m21

!print*,'loop5'
call write_symm_xmi_1(i7,m21,m5,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 705
if (sf2.eq.1)goto 100

!i5up6=i5
i6i7=i7
m6m21=m21

do m6=m5+1,nqul
!i5=i5up6
i7=i6i7
m21=m6m21

!print*,'loop6'
call write_symm_xmi_1(i7,m21,m6,sf1,sf2,str3,q_fac2)
!print*,'sf1,sf2',sf1,sf2
if (sf1.eq.1)goto 706
if (sf2.eq.1)goto 100

!i5up7=i5
i7i7=i7
m7m21=m21

do m7=m6+1,nqul
!i5=i5up7
i7=i7i7
m21=m7m21

!print*,'loop7'
call write_symm_xmi_1(i7,m21,m7,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 707
if (sf2.eq.1)goto 100

!i5up8=i5
i8i7=i7
m8m21=m21

do m8=m7+1,nqul
!i5=i5up8
i7=i8i7
m21=m8m21

!print*,'loop8'
call write_symm_xmi_1(i7,m21,m8,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 708
if (sf2.eq.1)goto 100

!i5up9=i5
i9i7=i7
m9m21=m21

do m9=m8+1,nqul
!i5=i5up9
i7=i9i7
m21=m9m21

!print*,'loop9'
call write_symm_xmi_1(i7,m21,m9,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 709
if (sf2.eq.1)goto 100

!i5up10=i5
i10i7=i7
m10m21=m21

do m10=m9+1,nqul
!i5=i5up10
i7=i10i7
m21=m10m21

!print*,'loop10'
call write_symm_xmi_1(i7,m21,m10,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 710
if (sf2.eq.1)goto 100

!i5up11=i5
i11i7=i7
m11m21=m21

do m11=m10+1,nqul
!i5=i5up11
i7=i11i7
m21=m11m21

!print*,'loop11'
call write_symm_xmi_1(i7,m21,m11,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 711
if (sf2.eq.1)goto 100

!i5up12=i5
i12i7=i7
m12m21=m21

do m12=m11+1,nqul
!i5=i5up12
i7=i12i7
m21=m12m21

print*,'loop12'
call write_symm_xmi_1(i7,m21,m12,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 712
if (sf2.eq.1)goto 100

!i5up13=i5
i13i7=i7
m13m21=m21

do m13=m12+1,nqul
!i5=i5up13
i7=i13i7
m21=m13m21

print*,'loop13'
call write_symm_xmi_1(i7,m21,m13,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 713
if (sf2.eq.1)goto 100

!i5up14=i5
i14i7=i7
m14m21=m21

do m14=m13+1,nqul
!i5=i5up14
i7=i14i7
m21=m14m21

!print*,'loop14',m14,sf1,sf2
!stop
call write_symm_xmi_1(i7,m21,m14,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 714
if (sf2.eq.1)goto 100

!i5up15=i5
i15i7=i7
m15m21=m21

do m15=m14+1,nqul
!i5=i5up15
i7=i15i7
m21=m15m21

print*,'loop15'
call write_symm_xmi_1(i7,m21,m15,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 715
if (sf2.eq.1)goto 100

!i5up16=i5
i16i7=i7
m16m21=m21

do m16=m15+1,nqul
!i5=i5up16
i7=i16i7
m21=m16m21

print*,'loop16'
call write_symm_xmi_1(i7,m21,m16,sf1,sf2,str3,q_fac2)
if (sf1.eq.1)goto 716
if (sf2.eq.1)goto 100

!i5up16=i5
i16i7=i7
m16m21=m21

716 enddo

715 enddo

714 enddo

713 enddo

712 enddo

711 enddo

710 enddo

709 enddo

708 enddo

707 enddo

706 enddo

705 enddo

704 enddo

703 enddo

702 enddo

701 enddo


close(21)

if(incmplt.eq.1)then
write(10,*)'----------- incomplete set -----------'
do i=1,mincmplt
 qq(i)=q_fac2(mincmplt_set(i))
 qq1(i)=str_quality_1(mincmplt_set(i))
 qq2(i)=str_quality_2(mincmplt_set(i))
 bondq4(i)=bondq(mincmplt_set(i))
    if(niao.eq.0)then
     write(10,900),qq1(i),bondq4(i),qq2(i),qq(i),'|',(str3(mincmplt_set(i),m20),m20=1,nae)
    endif
    if(niao.gt.1)then
     write(10,901),qq1(i),bondq4(i),qq2(i),qq(i),'|',1,':',niao,(str3(mincmplt_set(i),m20),m20=1,nae)
    endif
    if(niao.eq.1)then
     write(10,909),qq1(i),bondq4(i),qq2(i),qq(i),'|',1,1,(str3(mincmplt_set(i),m20),m20=1,nae)
    endif
enddo
endif

900 format(I3,x,I3,x,I3,x,I3,x,a,x,25I4)
901 format(I3,x,I3,x,I3,x,I3,x,a,x,I1,a,I3,x,25I4)
909 format(I3,x,I3,x,I3,x,I3,x,a,x,I3,I3,x,25I4)

open(unit=121,file='script10',status='unknown')
write(121,111)'rm -rf','Rumer_Sets.dat'
close(121)
CALL SYSTEM ("chmod +x script10 ")
CALL SYSTEM ("./script10 ")
CALL SYSTEM ("rm script10 ")
111 format(a,x,a)
100 return
end subroutine write_symm_xmi_new
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_symm_xmi_1(i7,m21,m3,sf1,sf2,str3,q_fac2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none

integer::nl,strn,ncqs,tostr,initstr,i,i1,i2,i3,i4,i5,i6,i7,i8,i9,m119,m18,m19,m20,m21,m23,m24,count&
,qul(100),nqul,j,jj,jjj,fg,flg,ii5,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,x,y,strs
integer::i5up16,i16i7,m16m21,i5up15,i15i7,m15m21,i5up14,i14i7,m14m21,i5up13,i13i7,m13m21,i5up12,i12i7,m12m21,&
i5up11,i11i7,m11m21,i5up10,i10i7,m10m21,i5up9,i9i7,m9m21,i5up8,i8i7,m8m21,i5up7,i7i7,m7m21,i5up6,i6i7,m6m21,&
i5up5,i5i7,m5m21,i5up4,i4i7,m4m21,i5up3,i3i7,m3m21,i5up2,i2i7,m2m21,i5up1,i1i7,m1m21,i5up17,i17i7,m17m21
integer::str3(15000,20),q_fac2(15000),finalvec(15000),strset(1000),col(1000),sigsym(15000),tnqs_sig,&
ffvec2(15000,1000),bondq(15000),bondq4(15000),nqset(15000),str5(2000,20),nstr7,&
tndet,totstr,Ifail,indpnt,strno(1000),str_quality_1(15000),str_quality_2(15000),ttqlty0,ttqlty&
,tqlty,bqlty,sqlty,hqlty,tnqs,nssym,qulsym(15000),symq(15000),set_number,ttqlty1,det_inv,ttqlty2,ttqlty3
integer::rumer(15000),rumer_rad(15000),quality_fac(15000),rumset,u1,max_set,mns,mm,group_num(15000)
!integer,dimension(:),allocatable::qq1,qq2,qq
integer::qq1(5000),qq2(5000),qq(5000),str2(2000,20),Rid,set_num(100),sf1,sf2,incmplt,mincmplt,mincmplt_set(500)
integer::j1,j2,grp(30)
real*8::ovlp
Double Precision::D(1000)
character(10)::dd,a
character(len=100)::outfile

common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig
common/str/str5,nstr7
common/infosymm/set_number,hqlty,ttqlty0,ttqlty1,ttqlty2,ttqlty3,ncqs,qul,Rid&
,mns,u1,max_set,rumset,nqset,strset,strn,totstr,strno,incmplt,mincmplt,mincmplt_set,group_num

print*,'*****************************'
print*,'enter write_symm_xmi_1'
print*,''
print*,''


!do i=i4+1,i5
!do i=i5+1,15000
!finalvec(i)=0
!enddo
!do m19=1,20
!print*,(str2(m19,i8),i8=1,8)
!enddo

!print*,'i7,m21',i7,m21
!stop
do m19=i7+1,1000
strno(m19)=0
 do i8=1,15
 str2(m19,i8)=0
 enddo
enddo

!i7=i3
totstr=i7
!m21=mm

jj=0
sf1=0
sf2=0
flg=0

!print*,'nqset,strn',nqset(m3),strn
if(nqset(m3).gt.strn)then
sf1=1
!print*,'rejected'
write(*,558)'number of structures of this group is =',nqset(m3),'so the group is rejected'
goto 100
endif
558 format(a,x,I2,x,a)

do i=strset(m3-1)+1,strset(m3)
totstr=totstr+1
strno(totstr)=i
!print*,'strno',strno(totstr)
enddo

j=0
do i=1,totstr
if(i.ne.1)then
j1=group_num(strno(i-1))
j2=group_num(strno(i))
if(j1.eq.j2)goto 555
endif
j=j+1
grp(j)=group_num(strno(i))
555 enddo
write(*,556),'groups involved = ',(grp(i),i=1,j)
556 format(a,3x,20I3)
write(*,*)

if(totstr.gt.strn)then
sf1=1
!print*,'rejected'
write(*,557)'includig this group the total number of structures &
becomes',totstr,'which is biggere than',strn,'so group has been rejected'
goto 100
endif
557 format(a,x,I2,x,a,x,I2,x,a)

write(*,559)'grps','str_num','structures'
do i=1,totstr
write(*,191)group_num(strno(i)),strno(i),(str3(strno(i),j),j=1,nae)
enddo
191 format(I3,3x,3x,I3,3X,20I3)
559 format(a,5x,a,6x,a)

call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)

if(Ifail.eq.0)print*,'structures are INDEPENDENT'
if(Ifail.eq.1)print*,'structures are DEPENDENT'

if(Ifail.eq.1) then
sf1=1
goto 100
endif

if(totstr.gt.mincmplt)then
mincmplt=totstr
do i=1,mincmplt
mincmplt_set(i)=strno(i)
enddo
endif

if(totstr.gt.1.and.ovopt.eq.vpt) then
 call MatLDR('str',strno,totstr,D)

 ovlp=1.0
 do i=1,totstr
 ovlp=ovlp*D(i)
 enddo

 if(1.0-ovlp.gt.ovval)then 
 sf1=1
 goto 100
 endif
endif

 m21=m21+nqset(m3)
 if(m21.gt.strn)then
 sf1=1
 goto 100
 endif

 do m19=strset(m3-1)+1,strset(m3)
 if(q_fac2(m19).ne.qul(m3))goto 603
 i7=i7+1
! print*,'new struc',i7,m19,q_fac2(m19),str_quality_1(m19),str_quality_2(m19),bondq(m19)
  do i3=1,nae
  str2(i7,i3)=str3(m19,i3)
  enddo
 col(i7)=m19
 qq(i7)=q_fac2(m19)
 qq1(i7)=str_quality_1(m19)
 qq2(i7)=str_quality_2(m19)
 bondq4(i7)=bondq(m19)
 flg=1
 count=3
 603 enddo

 if(i7.eq.strn) then
  ttqlty=0
  tqlty=0
  bqlty=0
  sqlty=0
  do m119=1,i7
   ttqlty=ttqlty+qq(m119)
  enddo
!   print*,'ttqlty,ttqlty0',ttqlty,ttqlty0
  if(ttqlty.le.ttqlty0)then
    set_number=set_number+1
    mns=mns+1
    if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
    if(nfset.eq.1.and.set_number.gt.1)then
     if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!     write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
    endif

print*,'i7***********',i7,m21,nfset
   do m19=1,i7
incmplt=0
    if(niao.eq.0)then
     write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
    endif
    if(niao.gt.1)then
     write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
    endif
    if(niao.eq.1)then
     write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
    endif
   enddo

   if(ovopt.eq.vpt) then
    call MatLDR('str',col,i7,D)

    ovlp=1.0
    do i=1,i7
     ovlp=ovlp*D(i)
    enddo
    write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
   endif

    if(nfset.eq.3)then
    Rumwrite=1
    call Rumer_set_id(str2,i7,nl,Rid,set_num)
    endif

   if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
   if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
   bqlty=0
   write(9+u1,*),'Set_number=',set_number
   write(9+u1,*),'    '
  endif

  if(nfset.eq.0.and.ttqlty.le.ttqlty0)then 
   sf2=1
   sf1=0
   goto 100
  endif
! endif

print*,'i am here',mns,max_set,u1,mset,mset-1
 if(mns.eq.max_set)then
  if(u1.eq.mset-1)then
print*,'i am here',mns,max_set,u1,mset,mset-1
   sf2=0
   sf1=1
   goto 100
  endif
  close(9+u1)
  mns=0
  u1=u1+1
  if(u1.le.9)then
   write(a,'(I1)')u1
  endif
  if(u1.gt.9.and.u1.lt.100)then
   write(a,'(I2)')u1
  endif
  if(u1.gt.99.and.u1.lt.1000)then
   write(a,'(I3)')u1
  endif
  if(u1.eq.1000)then
   write(*,*)'Maximum output file has been set to 1000, it seems that you have more&
   sets. Please increase the limit'
   sf2=1
   goto 100
  endif
  outfile='structure_set_'//trim(a)//trim('.dat')
  open(unit=9+u1,file=outfile,status='unknown')
 endif

!sf2=1
goto 100
endif

900 format(I3,x,I3,x,I3,x,I3,x,a,x,25I4)
901 format(I3,x,I3,x,I3,x,I3,x,a,x,I1,a,I3,x,25I4)
909 format(I3,x,I3,x,I3,x,I3,x,a,x,I3,I3,x,25I4)
!910 format(a,x,a,x,a,x,a)
!911 format(15x,I3,7x,I3,7x,I3)
910 format(a,I3)
920 format(a,2x,I5,2x,a,2x,100I5)
911 format(15x,I3,7x,I3,7x,I3)
912 format(a,3x,F10.3)

100 return
end subroutine write_symm_xmi_1
