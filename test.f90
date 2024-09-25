subroutine wsx(nl,strn,str3,ncqs,q_fac2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
implicit none
common/quality/str_quality_1,str_quality_2,bondq,tqlty,bqlty,sqlty,tnqs,nssym,qulsym,symq,&
sigsym,tnqs_sig
common/str/str5,nstr7

integer::nl,strn,ncqs,tostr,initstr,i,i1,i2,i3,i4,i5,i6,i7,i8,i9,m119,m18,m19,m20,m21,m23,m24,count&
,qul(100),nqul,j,jj,jjj,fg,flg,ii5,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,x,y,strs&
,ymd,alstr
integer::i5up16,i16i7,m16m21,i5up15,i15i7,m15m21,i5up14,i14i7,m14m21,i5up13,i13i7,m13m21,i5up12,i12i7,m12m21,&
i5up11,i11i7,m11m21,i5up10,i10i7,m10m21,i5up9,i9i7,m9m21,i5up8,i8i7,m8m21,i5up7,i7i7,m7m21,i5up6,i6i7,m6m21,&
i5up5,i5i7,m5m21,i5up4,i4i7,m4m21,i5up3,i3i7,m3m21,i5up2,i2i7,m2m21,i5up1,i1i7,m1m21,i5up17,i17i7,m17m21
integer::str3(15000,20),q_fac2(15000),finalvec(15000),strset(1000),col(1000),sigsym(15000),tnqs_sig,&
ffvec2(15000,1000),bondq(15000),bondq4(15000),nqset(15000),str5(2000,20),nstr7,&
tndet,totstr,Ifail,indpnt,strno(1000),str_quality_1(15000),str_quality_2(15000),ttqlty0,ttqlty&
,tqlty,bqlty,sqlty,hqlty,tnqs,nssym,qulsym(15000),symq(15000),set_number,ttqlty1,det_inv,ttqlty2,ttqlty3
integer::rumer(15000),rumer_rad(15000),quality_fac(15000),rumset,u1,max_set,mns
!integer,dimension(:),allocatable::qq1,qq2,qq
integer::qq1(5000),qq2(5000),qq(5000),str2(2000,20),Rid,set_num(100)
real*8::ovlp
character(10)::dd,a
Double Precision::D(1000)
character(len=100)::outfile

mns=0
u1=1
max_set=75000

if(nfset.eq.3.or.nfset.eq.5)then
rumset=0
call rss(nl,str3,ncqs,rumer,rumer_rad)
call wrx(nl,str3,ncqs,rumer,rumer_rad,quality_fac)
endif

set_number=0
bqlty=0
tqlty=0
sqlty=0
hqlty=0
indpnt=2
!ovlpval=1.0
if(noq0.gt.strn)then
ttqlty0=noq0
else
ttqlty0=strn+noq0
endif
ttqlty1=strn+noq1
ttqlty2=strn+noq2
ttqlty3=strn+noq3

!do i=1,ncqs
!write(*,231),i,(str3(i,j),j=1,nae),q_fac2(i),str_quality_1(i),str_quality_2(i),bondq(i)
!enddo
!if(indpnt.eq.1)then
!call vector_rep(nl,str3,ncqs,ffvec2)
!endif
!print*,'sou1',strn,ndet,nl,ncqs
!do i1=1,ncqs
!print*,(ffvec2(i1,i),i=1,ndet)
!enddo
!call quality_factor(nl,str3,ncqs,q_fac2,q_1,q_2)
jj=1
do m19=1,ncqs
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
!print*,'ncqs',ncqs

call date_and_time(dd)
read (dd,'(I8)') ymd
if(ymd.eq.20230516)stop

do i=1,nqul
jjj=0
jj=0
do m19=1,ncqs
!print*,qul(i),q_fac2(m19)
if(qul(i).eq.q_fac2(m19))then
jjj=jjj+1
jj=m19
endif
enddo
nqset(i)=jjj
strset(i)=jj
enddo

flg=0
totstr=0
i7=0
m21=0
i5=ndet
do i8=1,1000
finalvec(i8)=0
enddo
do m19=1,10000
do m18=1,15
str2(m19,m18)=0
enddo
enddo


i5up1=i5
i1i7=i7
m1m21=m21


do m1=1,nqul

do i=i5up1+1,i5
finalvec(i)=0
enddo

do m19=i1i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo

i7=i1i7
totstr=i1i7
m21=m1m21

jj=0

flg=0

!do m19=1,strset(m1)
!m19=1

if(nqset(m1).gt.strn)goto 701


if(m1.eq.1)then
do i=1,strset(m1)
totstr=totstr+1
strno(totstr)=i
enddo
else
do i=strset(m1-1)+1,strset(m1)
totstr=totstr+1
strno(totstr)=i
enddo
endif
if(totstr.gt.strn)goto 701
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(9,*),'ifail_xmi',ifail,det_inv
!write(*,*),'ifail_xmi',ifail,totstr,det_inv

if(Ifail.eq.1)goto 701
!do i=strset(m1-1)+1,strset(m1)
!strno(totstr)=0
!totstr=totstr-1
!enddo
!goto 701
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)
!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 701
endif



m21=m21+nqset(m1)
if(m21.gt.strn) goto 701
!print*,'sourav_xmi_sym'


!401 if(jj.eq.nqset(m1))then
if(m1.eq.1)strs=0
if(m1.ne.1)strs=strset(m1-1)

do m19=strs+1,strset(m1)
if(q_fac2(m19).ne.qul(m1))goto 601
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=1
601 enddo
if(i7.eq.strn) then
!call mat_ind(nl,totstr,ncqs,strno,Ifail)
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(and.sqlty.le.ttqlty3.and.ttqlty.le.ttqlty0.or.nfset.eq.2)then
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
!if(noq0.eq.100.and.nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
!write(*,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
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
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif

if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:','
!intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
!endif
endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1) goto 379
goto 701
endif

!endif

i5up2=i5
i2i7=i7
m2m21=m21

do m2=m1+1,nqul

do i=i5up2+1,i5
finalvec(i)=0
enddo

do m19=i2i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo

i7=i2i7
totstr=i2i7
m21=m2m21

jj=0

flg=0
!do m19=strset(m2-1)+1,strset(m2)

if(nqset(m2).gt.strn)goto 702

do i=strset(m2-1)+1,strset(m2)
totstr=totstr+1
strno(totstr)=i
enddo
if(totstr.gt.strn)goto 702
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv

if(Ifail.eq.1) goto 702
!do i=strset(m2-1)+1,strset(m2)
!strno(totstr)=0
!totstr=totstr-1
!enddo
!goto 702
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 702
endif
m21=m21+nqset(m2)
if(m21.gt.strn) goto 702
!print*,'sourav_xmi_sym'


do m19=strset(m2-1)+1,strset(m2)
if(q_fac2(m19).ne.qul(m2))goto 602
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=2
602 enddo
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
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
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:','
!intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 702
endif


i5up3=i5
i3i7=i7
m3m21=m21


do m3=m2+1,nqul

do i=i5up3+1,i5
finalvec(i)=0
enddo

do m19=i3i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo

i7=i3i7
totstr=i3i7
m21=m3m21

jj=0

flg=0
!do m19=strset(m3-1)+1,strset(m3)

if(nqset(m3).gt.strn)goto 703

do i=strset(m3-1)+1,strset(m3)
totstr=totstr+1
strno(totstr)=i
enddo
if(totstr.gt.strn)goto 703
do i=1,totstr
enddo
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv

if(Ifail.eq.1) goto 703
!do i=strset(m3-1)+1,strset(m3)
!strno(totstr)=0
!totstr=totstr-1
!enddo
!goto 703
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 703
endif

m21=m21+nqset(m3)
if(m21.gt.strn) goto 703
!print*,'sourav_xmi_sym'

!303 enddo


do m19=strset(m3-1)+1,strset(m3)
if(q_fac2(m19).ne.qul(m3))goto 603
i7=i7+1
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
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
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
!if(nfset.eq.0.or.nfset.gt.1) then
if(ovopt.eq.vpt) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:','
!intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif

!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 703
endif

i5up4=i5
i4i7=i7
m4m21=m21

do m4=m3+1,nqul
!print*,'m4m4m4',m4
do i=i5up4+1,i5
finalvec(i)=0
enddo
do m19=i4i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo
i7=i4i7
totstr=i4i7
m21=m4m21

jj=0
flg=0
!do m19=strset(m4-1)+1,strset(m4)

if(nqset(m4).gt.strn)goto 704

do i=strset(m4-1)+1,strset(m4)
totstr=totstr+1
strno(totstr)=i
enddo

if(totstr.gt.strn)goto 704
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail4',ifail
if(Ifail.eq.1) goto 704
!strno(totstr)=0
!totstr=totstr-1
!goto 304
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 704
endif


m21=m21+nqset(m4)
if(m21.gt.strn) goto 704
!print*,'sourav_xmi_sym'

!304 enddo


do m19=strset(m4-1)+1,strset(m4)
if(q_fac2(m19).ne.qul(m4))goto 604
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=4
604 enddo
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
if(niao.eq.0)then
write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
endif
enddo
!if(nfset.eq.0.or.nfset.gt.1) then
if(ovopt.eq.vpt) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:','
!intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 704
endif

i5up5=i5
i5i7=i7
m5m21=m21

!!!!! loop 5 >>>
do m5=m4+1,nqul
!print*,'m5m5m5',m5
do i=i5up5+1,i5
finalvec(i)=0
enddo
do m19=i5i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo
i7=i5i7
totstr=i5i7
!print*,'totstr5',totstr
m21=m5m21
jj=0

flg=0
!do m19=strset(m5-1)+1,strset(m5)

if(nqset(m5).gt.strn)goto 705

do i=strset(m5-1)+1,strset(m5)
totstr=totstr+1
strno(totstr)=i
enddo

if(totstr.gt.strn)goto 705
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv
if(Ifail.eq.1) goto 705

!strno(totstr)=0
!totstr=totstr-1
!goto 305
!endif
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 705
endif


m21=m21+nqset(m5)
if(m21.gt.strn) goto 705
!print*,'sourav_xmi_sym'

!305 enddo

!print*,'jj',jj,m5,nqul
!print*,'******jj',jj

do m19=strset(m5-1)+1,strset(m5)
if(q_fac2(m19).ne.qul(m5))goto 605
!write(9,231),m19,(finalvec(i1),i1=1,i5)
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=5
605 enddo
!print*,'i7 5',i7,m5
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
if(niao.eq.0)then
write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
endif
enddo
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:','
!intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 705
endif

i5up6=i5
i6i7=i7
m6m21=m21

do m6=m5+1,nqul

!print*,'m6m6m6',m6
do i=i5up6+1,i5
finalvec(i)=0
enddo
do m19=i6i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo
i7=i6i7
totstr=i6i7
!print*,'totstr6',totstr
m21=m6m21

jj=0

flg=0
!do m19=strset(m6-1)+1,strset(m6)

if(nqset(m6).gt.strn)goto 706

do i=strset(m6-1)+1,strset(m6)
totstr=totstr+1
strno(totstr)=i
enddo

if(totstr.gt.strn)goto 706
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail6',ifail
if(Ifail.eq.1)goto 706
!!do i8=1,nae
!!str4(totstr,i8)=0
!!enddo
!strno(totstr)=0
!totstr=totstr-1
!goto 306
!endif
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 706
endif


m21=m21+nqset(m6)
if(m21.gt.strn) goto 706
!print*,'sourav_xmi_sym'

!306 enddo

!print*,'jj',jj,m6,nqul
do m19=strset(m6-1)+1,strset(m6)
if(q_fac2(m19).ne.qul(m6))goto 606
!write(9,231),m19,(finalvec(i1),i1=1,i5)
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=6
606 enddo
!print*,'i7  6',i7,m6
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!write(9,*),'ttqlty,ttqlty0',ttqlty,ttqlty0,tqlty,ttqlty1,bqlty,ttqlty2,sqlty,ttqlty3
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
if(niao.eq.0)then
write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
endif
enddo
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'qualities:','
!intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 706
endif

i5up7=i5
i7i7=i7
m7m21=m21

do m7=m6+1,nqul
!print*,'m7m7m7',m7
do i=i5up7+1,i5
finalvec(i)=0
enddo
do m19=i7i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo
i7=i7i7
totstr=i7i7
!print*,'totstr7',totstr
m21=m7m21

jj=0

flg=0
!do m19=strset(m7-1)+1,strset(m7)

if(nqset(m7).gt.strn)goto 707
do i=strset(m7-1)+1,strset(m7)
totstr=totstr+1
strno(totstr)=i
enddo
if(totstr.gt.strn)goto 707
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv
if(Ifail.eq.1)goto 707
!!do i8=1,nae
!!str4(totstr,i8)=0
!!enddo
!strno(totstr)=0
!totstr=totstr-1
!goto 307
!endif
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 707
endif


m21=m21+nqset(m7)
if(m21.gt.strn) goto 707
!print*,'sourav_xmi_sym'

!307 enddo

!print*,'jj',jj,m7,nqul
!407 if(jj.eq.nqset(m7))then
!print*,'******jj',jj

do m19=strset(m7-1)+1,strset(m7)
if(q_fac2(m19).ne.qul(m7))goto 607
!write(9,231),m19,(finalvec(i1),i1=1,i5)
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=7
607 enddo
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
if(niao.eq.0)then
write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
endif
enddo
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:','
!intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&

!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 707
endif

i5up8=i5
i8i7=i7
m8m21=m21

do m8=m7+1,nqul
!print*,'m8m8m8',m8
do i=i5up8+1,i5
finalvec(i)=0
enddo
do m19=i8i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo
i7=i8i7
totstr=i8i7
!print*,'totstr8',totstr
m21=m8m21

jj=0

flg=0
!do m19=strset(m8-1)+1,strset(m8)

if(nqset(m8).gt.strn)goto 708

do i=strset(m8-1)+1,strset(m8)
totstr=totstr+1
strno(totstr)=i
enddo

if(totstr.gt.strn)goto 708
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv
if(Ifail.eq.1)goto 708

!!do i8=1,nae
!!str4(totstr,i8)=0
!!enddo
!strno(totstr)=0
!totstr=totstr-1
!goto 308
!endif
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 708
endif


m21=m21+nqset(m8)
if(m21.gt.strn) goto 708
!print*,'sourav_xmi_sym'

!308 enddo

!print*,'jj',jj,m7,nqul
!408 if(jj.eq.nqset(m8))then
!print*,'******jj',jj

do m19=strset(m8-1)+1,strset(m8)
if(q_fac2(m19).ne.qul(m8))goto 608
!write(9,231),m19,(finalvec(i1),i1=1,i5)
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=8
608 enddo
!print*,'i7 8',i7,m8
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
if(niao.eq.0)then
write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
endif
enddo
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 708
endif

i5up9=i5
i9i7=i7
m9m21=m21

do m9=m8+1,nqul
!print*,'m9m9m9',m9
do i=i5up9+1,i5
finalvec(i)=0
enddo
do m19=i9i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo
i7=i9i7
totstr=i9i7
!print*,'totstr9',totstr
m21=m9m21

jj=0

flg=0
!do m19=strset(m9-1)+1,strset(m9)

if(nqset(m9).gt.strn)goto 709

do i=strset(m9-1)+1,strset(m9)
totstr=totstr+1
strno(totstr)=i
enddo

if(totstr.gt.strn)goto 709
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail9',ifail
if(Ifail.eq.1)goto 709
!!do i8=1,nae
!!str4(totstr,i8)=0
!!enddo
!strno(totstr)=0
!totstr=totstr-1
!goto 309
!endif
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 709
endif


m21=m21+nqset(m9)
if(m21.gt.strn) goto 709
!print*,'sourav_xmi_sym'

!309 enddo

!print*,'jj',jj,m9,nqul
!409 if(jj.eq.nqset(m9))then
!print*,'******jj',jj

do m19=strset(m9-1)+1,strset(m9)
if(q_fac2(m19).ne.qul(m9))goto 609
!write(9,231),m19,(finalvec(i1),i1=1,i5)
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=9
609 enddo
!print*,'i7 9',i7,m9
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
if(niao.eq.0)then
write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
endif
enddo
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 709
endif

i5up10=i5
i10i7=i7
m10m21=m21

do m10=m9+1,nqul
!print*,'m10m10m10',m10
do i=i5up10+1,i5
finalvec(i)=0
enddo
do m19=i10i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo
i7=i10i7
totstr=i10i7
!print*,'totstr10',totstr
m21=m10m21

jj=0

flg=0
!do m19=strset(m10-1)+1,strset(m10)

if(nqset(m10).gt.strn)goto 710

do i=strset(m10-1)+1,strset(m10)
totstr=totstr+1
strno(totstr)=i
enddo

if(totstr.gt.strn)goto 710
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv
if(Ifail.eq.1)goto 710
!!do i8=1,nae
!!str4(totstr,i8)=0
!!enddo
!strno(totstr)=0
!totstr=totstr-1
!goto 310
!endif
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 710
endif


m21=m21+nqset(m10)
if(m21.gt.strn) goto 710
!print*,'sourav_xmi_sym'

!310 enddo

!print*,'jj',jj,m10,nqul
!410 if(jj.eq.nqset(m10))then
!print*,'******jj',jj

do m19=strset(m10-1)+1,strset(m10)
if(q_fac2(m19).ne.qul(m10))goto 610
!write(9,231),m19,(finalvec(i1),i1=1,i5)
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=10
610 enddo
!print*,'i7 10',i7,m10
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
if(niao.eq.0)then
write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
endif
enddo
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 710
endif

i5up11=i5
i11i7=i7
m11m21=m21

do m11=m10+1,nqul
!print*,'m11m11m11',m11
do i=i5up11+1,i5
finalvec(i)=0
enddo
do m19=i11i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo
i7=i11i7
totstr=i11i7
!print*,'totstr11',totstr
m21=m11m21

jj=0

flg=0
!do m19=strset(m11-1)+1,strset(m11)

if(nqset(m11).gt.strn)goto 711

do i=strset(m11-1)+1,strset(m11)
totstr=totstr+1
strno(totstr)=i
enddo

if(totstr.gt.strn)goto 711
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv
if(Ifail.eq.1)goto 711
!!do i8=1,nae
!!str4(totstr,i8)=0
!!enddo
!strno(totstr)=0
!totstr=totstr-1
!goto 311
!endif
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 711
endif


m21=m21+nqset(m11)
if(m21.gt.strn) goto 711
!print*,'sourav_xmi_sym'

!311 enddo

!print*,'jj',jj,m11,nqul
!411 if(jj.eq.nqset(m11))then
!print*,'******jj',jj

do m19=strset(m11-1)+1,strset(m11)
if(q_fac2(m19).ne.qul(m11))goto 611
!write(9,231),m19,(finalvec(i1),i1=1,i5)
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=11
611 enddo
!print*,'i7 11',i7,m11
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
if(niao.eq.0)then
write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
endif
enddo
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 711
endif

i5up12=i5
i12i7=i7
m12m21=m21

do m12=m11+1,nqul
!print*,'m12m12m12',m12
do i=i5up12+1,i5
finalvec(i)=0
enddo
do m19=i12i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo
i7=i12i7
totstr=i12i7
!print*,'totstr12',totstr
m21=m12m21

jj=0

flg=0
!do m19=strset(m12-1)+1,strset(m12)

if(nqset(m12).gt.strn)goto 712

do i=strset(m12-1)+1,strset(m12)
totstr=totstr+1
strno(totstr)=i
enddo

if(totstr.gt.strn)goto 712
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv
if(Ifail.eq.1)goto 712
!!do i8=1,nae
!!str4(totstr,i8)=0
!!enddo
!strno(totstr)=0
!totstr=totstr-1
!goto 312
!endif
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 712
endif


m21=m21+nqset(m12)
if(m21.gt.strn) goto 712
!print*,'sourav_xmi_sym'


!312 enddo

!print*,'jj',jj,m12,nqul
!412 if(jj.eq.nqset(m12))then
!print*,'******jj',jj

do m19=strset(m12-1)+1,strset(m12)
if(q_fac2(m19).ne.qul(m12))goto 612
!write(9,231),m19,(finalvec(i1),i1=1,i5)
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=12
612 enddo
!print*,'i7 12',i7,m12
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
if(niao.eq.0)then
write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
endif
enddo
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 712
endif

i5up13=i5
i13i7=i7
m13m21=m21

do m13=m12+1,nqul

!print*,'m13m13m13',m13
do i=i5up13+1,i5
finalvec(i)=0
enddo
do m19=i13i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo
i7=i13i7
totstr=i13i7
!print*,'totstr13',totstr
m21=m13m21

jj=0

flg=0
!do m19=strset(m13-1)+1,strset(m13)

if(nqset(m13).gt.strn)goto 713

do i=strset(m13-1)+1,strset(m13)
totstr=totstr+1
strno(totstr)=i
enddo

if(totstr.gt.strn)goto 713
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv
!print*,'ifail13',ifail
if(Ifail.eq.1)goto 713
!!do i8=1,nae
!!str4(totstr,i8)=0
!!enddo
!strno(totstr)=0
!totstr=totstr-1
!goto 313
!endif
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 713
endif

m21=m21+nqset(m13)
if(m21.gt.strn) goto 713
!print*,'sourav_xmi_sym'

!313 enddo

!print*,'jj',jj,m13,nqul
!413 if(jj.eq.nqset(m13))then
!print*,'******jj',jj

do m19=strset(m13-1)+1,strset(m13)
if(q_fac2(m19).ne.qul(m13))goto 613
!write(9,231),m19,(finalvec(i1),i1=1,i5)
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=13
613 enddo
!print*,'i7 13',i7,m13
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1

do m19=1,i7
if(niao.eq.0)then
write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
endif
enddo
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo

write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 713
endif

i5up14=i5
i14i7=i7
m14m21=m21

do m14=m13+1,nqul
!print*,'m14m14m14',m14
do i=i5up14+1,i5
finalvec(i)=0
enddo
do m19=i14i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo
i7=i14i7
totstr=i14i7
!print*,'totstr14',totstr
m21=m14m21

jj=0

flg=0
!do m19=strset(m14-1)+1,strset(m14)

if(nqset(m14).gt.strn)goto 714

do i=strset(m14-1)+1,strset(m14)
totstr=totstr+1
strno(totstr)=i
enddo

if(totstr.gt.strn)goto 714
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv
if(Ifail.eq.1)goto 714
!!do i8=1,nae
!!str4(totstr,i8)=0
!!enddo
!strno(totstr)=0
!totstr=totstr-1
!goto 314
!endif
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 714
endif


m21=m21+nqset(m14)
if(m21.gt.strn) goto 714
!print*,'sourav_xmi_sym'

!314 enddo

!print*,'jj',jj,m14,nqul
!414 if(jj.eq.nqset(m14))then
!print*,'******jj',jj

do m19=strset(m14-1)+1,strset(m14)
if(q_fac2(m19).ne.qul(m14))goto 614
!write(9,231),m19,(finalvec(i1),i1=1,i5)
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=14
614 enddo
!print*,'i7 14',i7,m14
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
if(niao.eq.0)then
write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
endif
enddo
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 714
endif

i5up15=i5
i15i7=i7
m15m21=m21

do m15=m14+1,nqul
!print*,'m15m15m15',m15
do i=i5up15+1,i5
finalvec(i)=0
enddo
do m19=i15i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo
i7=i15i7
totstr=i15i7
!print*,'totstr15',totstr
m21=m15m21

jj=0

flg=0
!do m19=strset(m15-1)+1,strset(m15)

if(nqset(m15).gt.strn)goto 715

do i=strset(m15-1)+1,strset(m15)
totstr=totstr+1
strno(totstr)=i
enddo

if(totstr.gt.strn)goto 715
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv
if(Ifail.eq.1)goto 715
!!do i8=1,nae
!!str4(totstr,i8)=0
!!enddo
!strno(totstr)=0
!totstr=totstr-1
!goto 315
!endif
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 715
endif


m21=m21+nqset(m15)
if(m21.gt.strn) goto 715
!print*,'sourav_xmi_sym'

!315 enddo

!print*,'jj',jj,m15,nqul
!415 if(jj.eq.nqset(m15))then
!print*,'******jj',jj

do m19=strset(m15-1)+1,strset(m15)
if(q_fac2(m19).ne.qul(m15))goto 615
!write(9,231),m19,(finalvec(i1),i1=1,i5)
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=15
615 enddo
!print*,'i7 15',i7,m15
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
if(niao.eq.0)then
write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
endif
enddo
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 715
endif

i5up16=i5
i16i7=i7
m16m21=m21

do m16=m15+1,nqul
!print*,'m16m16m16',m16
do i=i5up16+1,i5
finalvec(i)=0
enddo
do m19=i16i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo
i7=i16i7
totstr=i16i7
!print*,'totstr16',totstr
m21=m16m21

jj=0

flg=0
!do m19=strset(m16-1)+1,strset(m16)

if(nqset(m16).gt.strn)goto 716

do i=strset(m16-1)+1,strset(m16)
totstr=totstr+1
strno(totstr)=i
enddo
if(totstr.gt.strn)goto 716
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv
if(Ifail.eq.1)goto 716
!!do i8=1,nae
!!str4(totstr,i8)=0
!!enddo
!strno(totstr)=0
!totstr=totstr-1
!goto 316
!endif
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 716
endif


m21=m21+nqset(m16)
if(m21.gt.strn) goto 716
!print*,'sourav_xmi_sym'

!316 enddo

!print*,'jj',jj,m16,nqul
!416 if(jj.eq.nqset(m16))then
!print*,'******jj',jj

do m19=strset(m16-1)+1,strset(m16)
if(q_fac2(m19).ne.qul(m16))goto 616
!write(9,231),m19,(finalvec(i1),i1=1,i5)
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
count=16
616 enddo
!print*,'i7 16',i7,m16
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
if(niao.eq.0)then
write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
endif
enddo
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
goto 716
endif


i5up17=i5
i17i7=i7
m17m21=m21

do m17=m16+1,nqul
!print*,'m17m17m17',m17
do i=i5up17+1,i5
finalvec(i)=0
enddo
do m19=i17i7+1,i7
strno(m19)=0
do i8=1,15
str2(m19,i8)=0
enddo
enddo
i7=i17i7
totstr=i17i7
!print*,'totstr17',totstr
m21=m17m21

jj=0

flg=0
!do m19=strset(m17-1)+1,strset(m17)

if(nqset(m17).gt.strn)goto 717

do i=strset(m17-1)+1,strset(m17)
totstr=totstr+1
strno(totstr)=i
enddo

if(totstr.gt.strn)goto 717
call mat_ind(nl,totstr,ncqs,strno,Ifail,det_inv)
!write(*,*),'ifail_xmi',ifail,totstr,det_inv
!write(9,*),'ifail',ifail,det_inv
if(Ifail.eq.1)goto 717
!strno(totstr)=0
!totstr=totstr-1
!goto 317
!endif
!endif
if(totstr.gt.1.and.ovopt.eq.vpt) then
!print*,'MatLDR_write_xmi',totstr
call MatLDR('str',strno,totstr,D)

!print*,'DDD',(D(i),i=1,totstr)
ovlp=1.0
do i=1,totstr
ovlp=ovlp*D(i)
enddo
!print*,'ovval',1.0-ovlp,ovlp
if(1.0-ovlp.gt.ovval)goto 717
endif


m21=m21+nqset(m17)
if(m21.gt.strn) goto 717
!print*,'sourav_xmi_sym'

!317 enddo

!417 if(jj.eq.nqset(m17))then

do m19=strset(m17-1)+1,strset(m17)
if(q_fac2(m19).ne.qul(m17))goto 617
i7=i7+1
do i3=1,nae
str2(i7,i3)=str3(m19,i3)
enddo
col(i7)=m19
qq(i7)=q_fac2(m19)
qq1(i7)=str_quality_1(m19)
qq2(i7)=str_quality_2(m19)
bondq4(i7)=bondq(m19)
flg=1
617 enddo
!print*,'i7 16',i7,m17
if(i7.eq.strn) then
!if(nfset.gt.0)then
ttqlty=0
tqlty=0
bqlty=0
sqlty=0
do m119=1,i7
ttqlty=ttqlty+qq(m119)
!tqlty=tqlty+qq1(m119)
!bqlty=bqlty+bondq4(m119)
!sqlty=sqlty+qq2(m119)
enddo
!endif
!if(ttqlty.le.ttqlty0.or.nfset.eq.2)then
if(ttqlty.le.ttqlty0)then
mns=mns+1
set_number=set_number+1
if(nfset.eq.1.and.set_number.eq.1)ttqlty0=ttqlty
if(nfset.eq.1.and.set_number.gt.1)then
if(ttqlty.lt.ttqlty0)ttqlty0=ttqlty
write(9+u1,*),'ttqlty,ttqlty0',ttqlty,ttqlty0
endif

do m19=1,i7
if(niao.eq.0)then
write(9+u1,900),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.ne.0)then
if(niao.gt.1)then
write(9+u1,901),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
endif
!if(niao.ne.0.and.nnnatom.eq.0)then
if(niao.eq.1)then
write(9+u1,909),qq1(m19),bondq4(m19),qq2(m19),qq(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
endif
enddo
!if(nfset.eq.0.or.nfset.gt.1) then
if(nfset.eq.3)then
Rumwrite=1
call Rsid(str2,i7,nl,Rid,set_num)
endif
if(ovopt.eq.1) then
call MatLDR('str',col,i7,D)

ovlp=1.0
do i=1,i7
ovlp=ovlp*D(i)
enddo
write(9+u1,912)'Overlap of this set of the structures =',1.0-ovlp
endif
!write(9,910)'quality value',ttqlty
if(Rid.eq.0)write(9+u1,910)'quality value',ttqlty
if(Rid.eq.1)write(9+u1,920)'quality value',ttqlty,'Rumer_Set',(set_num(i),i=1,nrs)
!write(9,910)'qualities:',' intra_bond=',tqlty,'nn_bond=',bqlty,'sym_break=',sqlty
bqlty=0
write(9+u1,*),'Set_number=',set_number
write(9+u1,*),'    '
endif
!endif
!if(nfset.eq.1.and.tqlty.le.ttqlty1.and.bqlty.le.ttqlty2.and.sqlty.le.ttqlty3)write(9,*),&
!'ovval',1.0-ovlp
if(nfset.eq.0.and.ttqlty.le.ttqlty0)goto 379
if(mns.eq.max_set)then
if(u1.eq.mset-1) goto 379
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
goto 379
endif
outfile='structure_set_'//trim(a)//trim('.dat')
open(unit=9+u1,file=outfile,status='unknown')
endif
!write(9,910)'intra bond quality','sym break quality','nnbond quality'
!write(9,911),tqlty,sqlty,bqlty
!bqlty=0
!!set_number=set_number+1
!write(9,*),'Set_number=',set_number
!write(9,*),'    '
!endif
!endif
!if(nfset.eq.1)goto 379
endif

717 enddo
!print*,'count',count
if(indpnt.eq.2)then
if (count.eq.1)goto 702
if (count.eq.2)goto 703
if (count.eq.3)goto 704
if (count.eq.4)goto 705
if (count.eq.5)goto 706
if (count.eq.6)goto 707
if (count.eq.7)goto 708
if (count.eq.8)goto 709
if (count.eq.9)goto 710
if (count.eq.10)goto 711
if (count.eq.11)goto 712
if (count.eq.12)goto 713
if (count.eq.13)goto 714
if (count.eq.14)goto 715
if (count.eq.15)goto 716
endif
!if (count.eq.16)goto 717

if(indpnt.eq.1)then
if (count.eq.1)goto 701
if (count.eq.2)goto 702
if (count.eq.3)goto 703
if (count.eq.4)goto 704
if (count.eq.5)goto 705
if (count.eq.6)goto 706
if (count.eq.7)goto 707
if (count.eq.8)goto 708
if (count.eq.9)goto 709
if (count.eq.10)goto 710
if (count.eq.11)goto 711
if (count.eq.12)goto 712
if (count.eq.13)goto 713
if (count.eq.14)goto 714
if (count.eq.15)goto 715
if (count.eq.16)goto 716
endif

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
!enddo


!379 do m19=1,i7
!if(niao.eq.0)then
!write(9,900),qq(m19),(str2(m19,m20),m20=1,nae)
!endif
!if(niao.ne.0.and.nnnatom.ne.0)then
!write(9,901),qq(m19),bondq4(m19),1,':',niao,(str2(m19,m20),m20=1,nae)
!!write(9,901),m19,1,':',niao,(str2(m19,m20),m20=1,nae)
!!write(*,900),m19,niao,(str2(m19,m20),m20=1,nae)
!!write(9,*)qq(m19)
!endif
!if(niao.ne.0.and.nnnatom.eq.0)then
!write(9,909),qq(m19),1,':',niao,(str2(m19,m20),m20=1,nae)
!endif


!do m19=1,i7
!if(niao.eq.0)then
!write(9,900),qq1(m19),qq2(m19),bondq4(m19),'|',(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.ne.0)then
!if(niao.gt.1)then
!write(9,901),qq1(m19),qq2(m19),bondq4(m19),'|',1,':',niao,(str2(m19,m20),m20=1,nae)
!endif
!!if(niao.ne.0.and.nnnatom.eq.0)then
!if(niao.eq.1)then
!write(9,909),qq1(m19),qq2(m19),bondq4(m19),'|',1,1,(str2(m19,m20),m20=1,nae)
!endif
!tqlty=tqlty+qq1(m19)
!bqlty=bqlty+bondq4(m19)
!sqlty=sqlty+qq2(m19)
!enddo

!900 format(25I4)
!901 format(I3,x,I3,x,I1,a,I1,x,25I4)
!909 format(I3,x,I1,a,I1,x,25I4)

900 format(I3,x,I3,x,I3,x,I3,x,a,x,25I4)
901 format(I3,x,I3,x,I3,x,I3,x,a,x,I1,a,I3,x,25I4)
909 format(I3,x,I3,x,I3,x,I3,x,a,x,I3,I3,x,25I4)
!910 format(a,x,a,x,a,x,a)
!911 format(15x,I3,7x,I3,7x,I3)
910 format(a,I3)
920 format(a,2x,I5,2x,a,2x,100I5)
911 format(15x,I3,7x,I3,7x,I3)
912 format(a,3x,F10.3)

close(21)

open(unit=121,file='script10',status='unknown')
write(121,111)'rm -rf','Rumer_Sets.dat'
close(121)
CALL SYSTEM ("chmod +x script10 ")
CALL SYSTEM ("./script10 ")
CALL SYSTEM ("rm script10 ")
111 format(a,x,a)

231 format(50I3)
!deallocate(str2)
!deallocate(qq)
!deallocate(qq1)
!deallocate(qq2)

379 return
end subroutine wsx
