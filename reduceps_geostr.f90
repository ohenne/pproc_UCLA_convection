! Compile (blizzard):
! xlf95_r -o reduceps_geostr reduceps_geostr.f90 `/sw/aix61/netcdf-4.2.1.1/bin/nf-config --fflags --flibs` 
! Compile (thunder):
! ifort -o reduceps_geostr reduceps_geostr.f90 -L/sw/sles10-x64/netcdf-3.6.3-intel/lib -lnetcdf -I/sw/sles10-x64/netcdf-3.6.3-intel/include

        Program reduce

        use netcdf
        implicit none
! This program glues together the ps files from UCLA
! Compile with ifort -o reduceps reduceps.f90 -L/sw/sles9-x64/netcdf-3.6.2-intel/lib -lnetcdf
! -I/sw/sles9-x64/netcdf-3.6.2-intel/include on tornado
! Assume number of time steps nt smaller than 10000, number of levels nlev smaller than 150
! and number of variables nv smaller than 100 ! otherwise need to be changed
        integer, parameter :: nv=200,nma=7,nmi=2,nsu=27,nde=9
	      integer :: nt, nlev, zid

        character(100) stem,nm,nm2
        character(200) outpath
        character(20) pref,name,dimname,cnx,cny
        character(4) xsuf,ysuf
        character(1) dum1
        character(2) dum2
        integer nx,ny
        integer status,ncid,ndim,nvar,ncidw,timeid,varid,ztid,zmid,xtype,ndims,dimids(5)
        integer flag(nv),flagz(nv)
        character(len=10) maxnms(nma)
        character(len=10) minnms(nmi)
        character(len=10) sumnms(nsu)
        character(len=10) delnms(nde)  ! Christopher
        character(len=100) namout(nv),lname(nv),uname(nv)
        real fval(nv)
        real, allocatable, dimension(:) :: time(:),zt(:),zm(:),dn0(:),fsttm(:),lsttm(:),nsmp(:)
        real, dimension(:,:,:), allocatable ::  var,varout
        integer, dimension(:,:,:), allocatable :: cnt
        integer i,j,k,kk,l, nlevr, ntr
!0. Some preparations
!--------------------
!* maxnms contains the names where we take max, minnms where we take min, summns where we take 
!*sum; if nothing specified will take the mean; list can be extended if needed

        maxnms(1)="cfl       "
        maxnms(2)="maxdiv    "
        maxnms(3)="wmax      "
        maxnms(4)="lmax     "
        maxnms(5)="bflxmx    "
        maxnms(6)="bflxrmx   "
        maxnms(7)="precip_m  "
        minnms(1)="bflxmn    "
        minnms(2)="bflxrmn   "
        sumnms(1)="wr_cs1    "
        sumnms(2)="wr_cs2    "
        sumnms(3)="wv_cs1    "
        sumnms(4)="wv_cs2    "
        sumnms(5)="wt_cs1    "
        sumnms(6)="wt_cs2    "
        sumnms(7)="rt_cs1    "
        sumnms(8)="rt_cs2    "
        sumnms(9)="rl_cs1    "
        sumnms(10)="rl_cs2    "
        sumnms(11)="tv_cs1    "
        sumnms(12)="tv_cs2    "
        sumnms(13)="tl_cs1    "
        sumnms(14)="tl_cs2    "
        sumnms(15)="w_cs1     "
        sumnms(16)="w_cs2     "
        sumnms(17)="cnt_cs1   "
        sumnms(18)="cnt_cs2   "
        !sumnms(19)="Nr        "  ! Christopher: Nr must now be averaged
        sumnms(19)="nrcnt     "
        sumnms(20)="hst_d1    "
        sumnms(21)="hst_d2    "
        sumnms(22)="hst_d3    "
        sumnms(23)="hst_d4    "
        sumnms(24)="dissum    "
        sumnms(25)="discnt    "
        sumnms(26)="dissum2   "
        sumnms(27)="discnt2   "
        delnms(1)="u_col      "  ! Christopher
        delnms(2)="v_col      "
        delnms(3)="w_col      "
        delnms(4)="t_col      "
        delnms(5)="q_col      "
        delnms(6)="l_col      "
        delnms(7)="r_col      "
        delnms(8)="n_col      "
        delnms(9)="p_col      "

!* Ask to get filenames
        call get_command_argument(1,stem)
        call get_command_argument(2,cnx)
        call get_command_argument(3,cny)
        call get_command_argument(4,outpath)
        read (cnx,*) nx
        read (cny,*) ny
!         print*,'Directory where files are'
!         read*,dirin
         print*,'Files prefix', stem
!         read*,pref
         print*,'Processors in x',nx
!         read*,nx
         print*,'Processors in y',ny
!         read*,ny
! Open the first file to check the amoutn of timesteps; allocate the arrays
        nm=trim(stem)//".ps.00000000.nc"
        print*,trim(nm)
        status=nf90_open(nm,nf90_nowrite,ncid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
        status=nf90_inquire(ncid,unlimitedDimID=timeid)
        status=nf90_inquire_dimension(ncid,timeid,len=nt)
        status=nf90_inq_dimid(ncid,"zt",zid)
        status=nf90_inquire_dimension(ncid,zid,len=nlev)
        status=nf90_close(ncid)
        print *, 'Number of timesteps ', nt
        allocate(var(nv,nlev,nt), varout(nv,nlev,nt), cnt(nv,nlev,nt))
        allocate(time(nt),zt(nlev),zm(nlev),dn0(nlev),fsttm(nt),lsttm(nt),nsmp(nt))
        varout(:,:,:)=0.
        cnt(:,:,:)=0
        var(:,:,:)=0.


!1. Open each file and read variable
!----------------------------------

!*Loop over files and reconstruct filename
        do j=1,ny
         if (j-1.lt.10) then
          write(dum1,'(I1.1)') j-1
          ysuf="000"//dum1
         else
          write(dum2,'(I2.1)') j-1
          ysuf="00"//dum2
         end if
         do i=1,nx
          if (i-1.lt.10) then
           write(dum1,'(I1.1)') i-1
           xsuf="000"//dum1
          else
           write(dum2,'(I2.1)') i-1
           xsuf="00"//dum2
          end if
          nm=trim(stem)//".ps."//ysuf//xsuf//".nc"
          print*,trim(nm)
           
!*Open file
          status=nf90_open(nm,nf90_nowrite,ncid)
           if (status.ne.nf90_noerr) print*,nf90_strerror(status)

!*Get Variable names, dims, attributes. Do this only once at it is the same for all files  
!*Also set the flag according to compute mean (0), max (2), min (3), sum (1)
          if ((i.eq.1).and.(j.eq.1)) then
           flag(:)=0
           status=nf90_inquire(ncid,ndim,nvar)
            if (status.ne.nf90_noerr) print*,nf90_strerror(status)
           status=nf90_inquire_dimension(ncid,1,dimname,ntr)
            if (status.ne.nf90_noerr) print*,nf90_strerror(status)
           status=nf90_inquire_dimension(ncid,2,dimname,nlevr)
            if (status.ne.nf90_noerr) print*,nf90_strerror(status)
           print*,'There are ',ntr,' time steps, ',nlevr,' levels and ',nvar, 'variables in the file'
           do k=1,nvar
            status=nf90_inquire_variable(ncid,k,name,xtype,ndims,dimids)
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
            namout(k)=name
!*Assume the attribute names are longname and units
            status=nf90_get_att(ncid,k,'longname',lname(k))
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
            status=nf90_get_att(ncid,k,'units',uname(k))
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
            status=nf90_get_att(ncid,k,'_FillValue',fval(k))
             if (status.ne.nf90_noerr) then 
               print*,nf90_strerror(status)
               fval(k) = -999.
             end if
!            print*,k,trim(namout(k)),' ',trim(lname(k)),' ',trim(uname(k))
            do kk=1,nma
              if (trim(name).eq.trim(maxnms(kk))) flag(k)=2
            end do
            do kk=1,nmi
              if (trim(name).eq.trim(minnms(kk))) flag(k)=3
            end do
            do kk=1,nsu
              if (trim(name).eq.trim(sumnms(kk))) flag(k)=1
            end do
            do kk=1,nde      ! Christopher
              if (trim(name).eq.trim(delnms(kk))) flag(k)=4
            end do

!*Need to remember whether the variable is on zt (flagz=0) or zm (flagz=1)
!*Assume the dimension order is time,zt and zm in file
            if (dimids(1).eq.2) flagz(k)=0
            if (dimids(1).eq.3) flagz(k)=1
           end do

!*Read time, zt, zm, dn0,fsttm,lsttm,nsmp  from first file only
!*Assume those are the only 1D-var and that they come in the following order in the file
            status=nf90_get_var(ncid,1,time(1:ntr))
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
            status=nf90_get_var(ncid,2,zt(1:nlevr))
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
            status=nf90_get_var(ncid,3,zm(1:nlevr))
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
            status=nf90_get_var(ncid,4,dn0(1:nlevr))
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
            status=nf90_get_var(ncid,5,fsttm(1:ntr))
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
            status=nf90_get_var(ncid,6,lsttm(1:ntr))
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
            status=nf90_get_var(ncid,7,nsmp(1:ntr))
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
          end if

!*Read remaining variable in var for all files
           !do k=10,nvar
           do k=8,nvar
	     if (flag(k).eq.4) cycle   ! Christopher
             status=nf90_get_var(ncid,k,var(k,1:nlevr,1:ntr))
             if (status.ne.nf90_noerr) print*,nf90_strerror(status)
           end do

!*Close file
          status=nf90_close(ncid)
           if (status.ne.nf90_noerr) print*,nf90_strerror(status)

!2. Compute max, min, sum or average
!----------------------------------
          !do k=10,nvar
          do k=8,nvar
           if (flag(k).le.1) then
            do kk=1,nlevr
            do l=1,ntr
             if (var(k,kk,l).gt.-999.) then
              varout(k,kk,l)=varout(k,kk,l)+var(k,kk,l)
              cnt(k,kk,l)=cnt(k,kk,l)+1
             end if
            end do
            end do
           else
            if (flag(k).eq.2) then
             do kk=1,nlevr
             do l=1,ntr
              if ((varout(k,kk,l).eq.0).or.(var(k,kk,l).gt.varout(k,kk,l))) varout(k,kk,l)=var(k,kk,l)
             end do
             end do
            else
             do kk=1,nlevr
             do l=1,ntr
              if ((varout(k,kk,l).eq.0).or.((var(k,kk,l).gt.-999).and.(var(k,kk,l).lt.varout(k,kk,l)))) varout(k,kk,l)=var(k,kk,l)
             end do
             end do
            endif
           end if
          end do

         end do
        end do
!*End of looping over the files
  
      do k=8,nvar
         if (flag(k).eq.0) varout(k,:,:)=varout(k,:,:)/cnt(k,:,:)
      end do

!3 Write new variables in *.ps.nc
!---------------------------------

!* Enter define mode
      nm2=trim(stem)//".ps.nc"
      print *,'a',nm2
      status=nf90_create(nm2, nf90_CLOBBER,ncidw)
       if (status.ne.nf90_noerr) print*,nf90_strerror(status)

!* Define dimenstions
      status=nf90_def_dim(ncidw,'time',nf90_unlimited,timeid)
       if (status.ne.nf90_noerr) print*,nf90_strerror(status)
      status=nf90_def_dim(ncidw,'zt',nlevr,ztid)
       if (status.ne.nf90_noerr) print*,nf90_strerror(status)
      status=nf90_def_dim(ncidw,'zm',nlevr,zmid)
       if (status.ne.nf90_noerr) print*,nf90_strerror(status)
!* Define and store variable
!* Care has to be taken for the first 7 (! Christopher) variables as they are one-dimensional
!* Assume the first 7 are time,zt,zm,dn0,fsstm,lsttm,nsmp
      do k=1,nvar
      if (flag(k).eq.4) cycle   ! Christopher
      ! Christopher: change units
      ! tot_qw (W/m^2 --> m/s g/kg)
      if (trim(namout(k)).eq.'tot_qw') then
        uname(k)='m/s g/kg'
        varout(k,1:nlevr,1:ntr) = varout(k,1:nlevr,1:ntr)*0.00032787
	print *,trim(namout(k))
      endif
      ! tot_lw (W/m^2 --> m/s g/kg)
      if (trim(namout(k)).eq.'tot_lw') then
        uname(k)='m/s g/kg'
        varout(k,1:nlevr,1:ntr) = varout(k,1:nlevr,1:ntr)*0.00032787
	print *,trim(namout(k))
      endif
      ! q_2 (kg^2/kg^2 --> g^2/kg^2)
      if (trim(namout(k)).eq.'q_2') then
        uname(k)='g^2/kg^2'
        varout(k,1:nlevr,1:ntr) = varout(k,1:nlevr,1:ntr)*1000000.
	print *,trim(namout(k))
      endif
      ! q_3 (kg^3/kg^3 --> g^3/kg^3)
      if (trim(namout(k)).eq.'q_3') then
        uname(k)='g^3/kg^3'
        varout(k,1:nlevr,1:ntr) = varout(k,1:nlevr,1:ntr)*1000000000.
	print *,trim(namout(k))
      endif
      ! l_2 (kg^2/kg^2 --> g^2/kg^2)
      if (trim(namout(k)).eq.'l_2') then
        uname(k)='g^2/kg^2'
        varout(k,1:nlevr,1:ntr) = varout(k,1:nlevr,1:ntr)*1000000.
	print *,trim(namout(k))
      endif
      ! l_3 (kg^3/kg^3 --> g^3/kg^3)
      if (trim(namout(k)).eq.'l_3') then
        uname(k)='g^3/kg^3'
        varout(k,1:nlevr,1:ntr) = varout(k,1:nlevr,1:ntr)*1000000000.
	print *,trim(namout(k))
      endif
      ! qt_th (kg/kg K --> g/kg K)
      if (trim(namout(k)).eq.'qt_th') then
        uname(k)='g/kg K'
        varout(k,1:nlevr,1:ntr) = varout(k,1:nlevr,1:ntr)*1000.
	print *,trim(namout(k))
      endif
      ! sfs_qw (W/m^2 --> m/s g/kg)
      if (trim(namout(k)).eq.'sfs_qw') then
        uname(k)='g/kg K'
        varout(k,1:nlevr,1:ntr) = varout(k,1:nlevr,1:ntr)*0.00032787
	print *,trim(namout(k))
      endif
      ! Christopher: end
      !if ((k.eq.1).or.(k.eq.7).or.(k.eq.8).or.(k.eq.9)) then
      if ((k.eq.1).or.(k.eq.5).or.(k.eq.6).or.(k.eq.7)) then
       status=nf90_def_var(ncidw,trim(namout(k)),NF90_float,(/timeid/),varid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
      end if
      if (k.eq.3) then
       status=nf90_def_var(ncidw,trim(namout(k)),NF90_float,(/zmid/),varid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
      end if
      !if ((k.eq.2).or.(k.eq.4).or.(k.eq.5).or.(k.eq.6)) then
      if ((k.eq.2).or.(k.eq.4)) then
       status=nf90_def_var(ncidw,trim(namout(k)),NF90_float,(/ztid/),varid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
      end if
      !if ((k.gt.9).and.(flagz(k).eq.0)) then
      if ((k.gt.7).and.(flagz(k).eq.0)) then
       status=nf90_def_var(ncidw,trim(namout(k)),NF90_float,(/ztid,timeid/),varid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
      end if
      !if ((k.gt.9).and.(flagz(k).eq.1)) then
      if ((k.gt.7).and.(flagz(k).eq.1)) then
       status=nf90_def_var(ncidw,trim(namout(k)),NF90_float,(/zmid,timeid/),varid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
      end if
       status=nf90_put_att(ncidw,varid,'longname',trim(lname(k)))
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
      status=nf90_put_att(ncidw,varid,'_FillValue',fval(k))
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
       status=nf90_put_att(ncidw,varid,'units',trim(uname(k)))
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
       status=nf90_enddef(ncid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
       !if (k.gt.9) then
       if (k.gt.7) then
        status=nf90_put_var(ncidw,varid,varout(k,1:nlevr,1:ntr))
         if (status.ne.nf90_noerr) print*,nf90_strerror(status)
       end if
       status=nf90_redef(ncid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
      end do
       status=nf90_enddef(ncid)

!*Store the 7(!) 1D fields
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)
        status=nf90_put_var(ncidw,1,time(1:ntr))
         if (status.ne.nf90_noerr) print*,nf90_strerror(status)
        status=nf90_put_var(ncidw,2,zt(1:nlevr))
         if (status.ne.nf90_noerr) print*,nf90_strerror(status)
        status=nf90_put_var(ncidw,3,zm(1:nlevr))
         if (status.ne.nf90_noerr) print*,nf90_strerror(status)
        status=nf90_put_var(ncidw,4,dn0(1:nlevr))
         if (status.ne.nf90_noerr) print*,nf90_strerror(status)
        status=nf90_put_var(ncidw,5,fsttm(1:ntr))
         if (status.ne.nf90_noerr) print*,nf90_strerror(status)
        status=nf90_put_var(ncidw,6,lsttm(1:ntr))
         if (status.ne.nf90_noerr) print*,nf90_strerror(status)
        status=nf90_put_var(ncidw,7,nsmp(1:ntr))
         if (status.ne.nf90_noerr) print*,nf90_strerror(status)

!* Close file
       status=nf90_close(ncid)
        if (status.ne.nf90_noerr) print*,nf90_strerror(status)

        end
