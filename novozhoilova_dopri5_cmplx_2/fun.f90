module fun
    use, intrinsic :: iso_c_binding
    use ifcore

    integer(c_int) ne, nt, nz, neqf, neqp, lworkf, lworkp, liworkf, liworkp, nrdf, nrdp, iparf, iparp, ioutf, ioutp, ididf, ididp, itolf, itolp
    real(c_double) zex, dz, tend, dtr(2), q(3), i(2), th(2), a(2), dcir(2), r(2), f0(3), dt, &
        pitch, f10, f20, f30, rtolf, rtolp, atolf, atolp, rparf, rparp, ftol, ptol, &
        x1r, x1i, q31, i1, r1, th1, dcir1, cos1, sin1, &
        x2r, x2i, q32, i2, r2, th2, dcir2, cos2, sin2, q3, a1, a2, ng
    complex(c_double_complex) fp(2)

    integer(c_int) breaknum(3)
    real(c_double) phitmp0(3), phitmp1(3)
    complex(c_double_complex) fc, fcomp(3)

    integer(c_int), allocatable, target :: idxre(:, :), idxim(:, :), iworkf(:), iworkp(:)
    complex(c_double_complex), allocatable, target :: mean(:), fcmplx(:, :)
    real(c_double), allocatable, target :: tax(:), zax(:), u(:), eta(:, :), etag(:, :), w(:, :), f(:, :), w1(:, :), p(:, :), &
                                           phi(:, :), phios(:, :), wos(:, :), workf(:), workp(:)

    complex(c_double_complex), parameter :: ic = (0.0d0, 1.0d0)
    real(c_double), parameter :: pi = 2.0D0*dacos(0.0D0)

    private freq_out, zex, tend, dtr, q, i, th, a, dcir, r, f0, pitch

contains

    subroutine init()
        implicit none

        integer(c_int) ii

        call read_param()

        nt = tend/dt + 1
        nz = zex/dz + 1

        neqp = 4*ne
        nrdp = 4*ne
        lworkp = 8*neqp + 5*nrdp + 21
        liworkp = nrdp + 21

        neqf = 6
        nrdf = 6
        lworkf = 11*neqf + 8*nrdf + 21
        liworkf = nrdf + 21

        call allocate_arrays()

        f(1, 1) = f10
        f(3, 1) = f20
        f(5, 1) = f30

        do ii = 1, nt
            tax(ii) = (ii - 1)*dt
        end do

        do ii = 1, nz
            zax(ii) = (ii - 1)*dz
        end do

        call calc_u(u, zex, nz, zax)

        do ii = 1, 2
            idxre(ii, :) = (/2*(ii - 1)*ne + 1:(2*ii - 1)*ne/)
            idxim(ii, :) = (/(2*ii - 1)*ne + 1:2*ii*ne/)
        end do
        
        ng = 2

        !for dfdtcmplx
        q31 = q(3)/q(1)
        i1 = i(1)
        r1 = r(1)
        th1 = th(1)
        dcir1 = dcir(1)

        q32 = q(3)/q(2)
        i2 = i(2)
        r2 = r(2)
        th2 = th(2)
        dcir2 = dcir(2)

        q3 = q(3)
        a1 = a(1)
        a2 = a(2)                

    end subroutine init

    subroutine allocate_arrays()
        use, intrinsic :: iso_c_binding
        implicit none

        integer(c_int) err_alloc

        allocate (f(6, nt), p(4*ne, nz), u(nz), tax(nt), zax(nz), mean(nz), eta(2, nt), etag(2, nt), w(3, nt - 1), w1(3, nt - 1), &
                  idxre(2, ne), idxim(2, ne), workf(lworkf), iworkf(liworkf), workp(lworkp), iworkp(liworkp), &
                  wos(3, nt - 1), phi(3, nt), phios(3, nt), fcmplx(3, nt), &
                  stat=err_alloc)

        if (err_alloc /= 0) then
            print *, "allocation error"
            pause
            stop
        end if
    end subroutine allocate_arrays

    subroutine deallocate_arrays()
        use, intrinsic :: iso_c_binding
        implicit none

        integer(c_int) err_dealloc

        deallocate (f, p, u, tax, zax, mean, eta, etag, w, w1, stat=err_dealloc)

        if (err_dealloc /= 0) then
            print *, "deallocation error"
            pause
            stop
        end if
    end subroutine deallocate_arrays

    subroutine read_param() bind(c, name='read_param')
        use, intrinsic :: iso_c_binding
        import
        implicit none

        namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, &
            dcir1, dcir2, r1, r2, f10, f20, f30, dt, dz, pitch, ftol, ptol

        real(c_double) q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, tol

        open (unit=1, file='input_fortran.dat', status='old', err=101)
        read (unit=1, nml=param, err=102)
        close (unit=1)

        q(1) = q1
        q(2) = q2
        q(3) = q3
        i(1) = i1
        i(2) = i2
        th(1) = th1
        th(2) = th2
        a(1) = a1
        a(2) = a2
        dtr(1) = dtr1
        dtr(2) = dtr2
        dcir(1) = dcir1
        dcir(2) = dcir2
        r(1) = r1
        r(2) = r2

        write (*, nml=param)

        return
101     print *, "error of file open"; pause; stop
102     print *, 'error of reading file "input_fortran.dat"'; pause; stop
    end subroutine read_param

    subroutine ode4f()
        import
        implicit none

        integer(c_int) i, j, aiparf(1), aiparp(1), itp, itf
        real(c_double) :: t, z, artolf(1), aatolf(1), arparf(1), artolp(1), aatolp(1), arparp(1), xoutp, xoutf
        real(c_double) yf(6), yp(neqp), pex(neqp)
        logical(4) pressed
        character(1) key
        integer(c_int), parameter :: esc = 27
        common/internp/xoutp, itp
        common/internf/xoutf, itf

        !solve eq. at t=0
        fp(1) = dcmplx(f(1,1), f(2,1))
        fp(2) = dcmplx(f(3,1), f(4,1))
       
        rparp = 0.0
        iparp = 0
        itolp = 0
        rtolp = ptol
        atolp = rtolp
        ioutp = neqp
        z = zax(1)
        xoutp = z
        itp = 0
        yp = p(:, 1)
        iworkp(:) = 0
        workp(:) = 0.0d0
        iworkp(5) = neqp

        artolp(1) = rtolp
        aatolp(1) = atolp
        arparp(1) = rparp
        aiparp(1) = iparp

        call dopri5_p(neqp, dpdz, z, yp, zex, artolp, aatolp, itolp, soloutp, ioutp, &
                      workp, lworkp, iworkp, liworkp, arparp, aiparp, ididp)

        p(:, nz) = yp(:)

        eta(:, 1) = eff(p(:, nz))
        etag(:, 1) = pitch**2/(pitch**2 + 1)*eta(:, 1)

        rparf = 0.0
        iparf = 0
        t = tax(1)
        xoutf = t
        itf = 0
        yf = f(:, 1)
        itolf = 0
        rtolf = ftol
        atolf = rtolf
        ioutf = 6

        artolf(1) = rtolf
        aatolf(1) = atolf
        arparf(1) = rparf
        aiparf(1) = iparf

        iworkf(:) = 0
        workf(:) = 0.0d0

        iworkf(5) = 6

        call dopri5_f(6, dfdt, t, yf, tend, artolf, aatolf, itolf, soloutf, ioutf, &
                      workf, lworkf, iworkf, liworkf, arparf, aiparf, ididf)

        do j = 1, neqf
            f(j, nt) = yf(j)
        end do
        call calcpex(f(:, nt), pex)
        eta(:, nt) = eff(pex)
        !eta(:, nt) = eff(p(:, nz))
        etag(:, nt) = pitch**2/(pitch**2 + 1)*eta(:, nt)

    end subroutine ode4f

    function eff(pex) result(eta)
        use, intrinsic :: iso_c_binding, only: c_double, c_int
        import, only:ne, idxre, idxim

        implicit none

        integer(c_int) i
        real(c_double) eta(2)
        real(c_double), intent(in) :: pex(:)

        do i = 1, 2
            eta(i) = 1 - sum(cdabs(dcmplx(pex(idxre(i, :)), pex(idxim(i, :))))**2)/ne
        end do
    end function eff

    subroutine dpdz(neqp, z, p, prhs, rparp, iparp)
        import :: ne, zex, f, ic, dtr
        implicit none

        real(c_double) z, p(*), prhs(*)

        integer(c_int) i, reidx(ne), imidx(ne), neqp, iparp
        real(c_double) u, rparp
        complex(c_double_complex) s(ne), ptmp(ne)

        u = dexp(-3*((z - zex/2)/(zex/2))**2)

        do i = 1, 2
            ptmp = dcmplx(p(idxre(i, :)), p(idxim(i, :)))

            s = ic*(fp(i)*u*dconjg(ptmp) - (dtr(i) + cdabs(ptmp)**2 - 1)*ptmp)

            prhs(idxre(i, :)) = dreal(s)
            prhs(idxim(i, :)) = dimag(s)
        end do
    end subroutine dpdz

    complex(c_double_complex) function xi(p, num)
        use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
        import, only:ne, nz, mean, u, dz, idxre, idxim

        implicit none

        integer(c_int) i, num
        real(c_double) p(:, :)

        do i = 1, nz
            mean(i) = sum(dcmplx(p(idxre(num, :), i), p(idxim(num, :), i))**2, 1)/ne
        end do

        mean = u*mean
        !mean = dconjg(u)*mean

        xi = (0.5d0*(mean(1) + mean(2)) + sum(mean(2:nz - 1)))*dz
    end function

    subroutine dfdt(neqf, t, f, s, rparf, iparf)
        implicit none

        integer(c_int) :: ii, jj, neqf, iparf, aiparp(1), itp
        real(c_double) t, z, f(neqf), yp(neqp), s(neqf), xoutp, &
            x1r, x1i, q31, i1, r1, th1, dcir1, cos1, sin1, &
            x2r, x2i, q32, i2, r2, th2, dcir2, cos2, sin2, q3, &
            f1, f2, f3, phi1, phi2, phi3, rparf, artolp(1), aatolp(1), arparp(1)
        complex(c_double_complex) x1, x2, fcmplx(neqf/2), scmplx(neqf/2)
        common/internp/xoutp, itp

        rparp = 0.0
        iparp = 0
        itolp = 0
        rtolp = ptol
        atolp = rtolp
        ioutp = neqp
        z = zax(1)
        xoutp = z
        itp = 0
        yp = p(:, 1)
        iworkp(:) = 0
        workp(:) = 0.0d0
        iworkp(5) = neqp

        artolp(1) = rtolp
        aatolp(1) = atolp
        arparp(1) = rparp
        aiparp(1) = iparp

        fp(1) = dcmplx(f(1), f(2))
        fp(2) = dcmplx(f(3), f(4))

        call dopri5_p(neqp, dpdz, z, yp, zex, artolp, aatolp, itolp, soloutp, ioutp, &
                      workp, lworkp, iworkp, liworkp, arparp, aiparp, ididp)

        p(:, nz) = yp(:)

        fcmplx(1) = dcmplx(f(1), f(2))
        fcmplx(2) = dcmplx(f(3), f(4))
        fcmplx(3) = dcmplx(f(5), f(6))

        scmplx = dfdtcmplx(fcmplx)

        s(1) = dreal(scmplx(1))
        s(2) = dimag(scmplx(1))
        s(3) = dreal(scmplx(2))
        s(4) = dimag(scmplx(2))
        s(5) = dreal(scmplx(3))
        s(6) = dimag(scmplx(3))
    end subroutine dfdt

    subroutine calc_u(u, zex, nz, zax)
        import
        implicit none

        integer(c_int), intent(in) :: nz
        real(c_double), intent(in) :: zex, zax(nz)
        real(c_double), intent(out) :: u(:)

        integer(c_int) i

        do i = 1, nz
            u(i) = dexp(-3*((zax(i) - zex/2)/(zex/2))**2)
        end do

    end subroutine

    subroutine soloutf(nr, xold, x, y, n, con, icomp, nd, rparf, iparf, irtrn)
        implicit none

        interface
            function contd5_f(ii, x, con, icomp, nd)
                implicit double precision(a - h, o - z)
                dimension con(5*nd), icomp(nd)
            end
        end interface

        integer(c_int) nr, n, nd, icomp(nd), iparf, irtrn, j, itf
        real(c_double) xold, x, con(5*nd), rparf, y(neqf), xoutf, pex(neqp), yy(neqf)
        complex(c_double_complex) fcmplx(3)
        logical(4) pressed
        character(1) key
        integer(c_int), parameter :: esc = 27
        common/internf/xoutf, itf

        if (nr .eq. 1) then
            itf = 1
            do j = 1, neqf
                f(j, itf) = y(j)
            end do
            call calcpex(y, pex)
            eta(:, itf) = eff(pex)
            !eta(:, itf) = eff(p(:, nz))
            etag(:, itf) = pitch**2/(pitch**2 + 1)*eta(:, itf)
            do j = 1, neqf/2
                fcmplx(j) = dcmplx(y(2*j - 1), y(2*j))
            end do
            write (*, '(a,f10.5,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,\,a)') 'Time = ', xoutf, '   |F1| = ', abs(fcmplx(1)), '   |F2| = ', abs(fcmplx(2)), &
                '   |F3| = ', abs(fcmplx(3)), '   Eff1 = ', eta(1, itf), '   Eff2 = ', eta(2, itf), char(13)
            xoutf = x + dt
        else
10          continue
            if (x .ge. xoutf) then
                itf = itf + 1
                do j = 1, neqf
                    f(j, itf) = contd5_f(j, xoutf, con, icomp, nd)
                end do
                call calcpex(f(:, itf), pex)
                eta(:, itf) = eff(pex)
                !eta(:, itf) = eff(p(:, nz))
                etag(:, itf) = pitch**2/(pitch**2 + 1)*eta(:, itf)
                do j = 1, neqf/2
                    fcmplx(j) = dcmplx(f(2*j - 1, itf), f(2*j, itf))
                end do
                write (*, '(a,f10.5,a,f10.6,a,f10.6,a,f10.6,a,f10.6,a,f10.6,\,a)') 'Time = ', xoutf, '   |F1| = ', abs(fcmplx(1)), '   |F2| = ', abs(fcmplx(2)), &
                    '   |F3| = ', abs(fcmplx(3)), '   Eff1 = ', eta(1, itf), '   Eff2 = ', eta(2, itf), char(13)
                xoutf = xoutf + dt
                goto 10
            end if
        end if

        pressed = peekcharqq()
        if (pressed) then
            key = getcharqq()
            if (ichar(key) .eq. esc) then
                write (*, '(/,a)') 'Quit?'
                key = getcharqq()
                if (ichar(key) .eq. 121 .or. ichar(key) .eq. 89) then
                    nt = itf
                    irtrn = -1
                    !return
                end if
            end if
        end if
        return
    end subroutine soloutf

    subroutine soloutp(nr, xold, x, y, n, con, icomp, nd, rparp, iparp, irtrn)
        implicit none

        interface
            function contd5_p(ii, x, con, icomp, nd)
                implicit double precision(a - h, o - z)
                dimension con(5*nd), icomp(nd)
            end
        end interface

        integer(c_int) nr, n, nd, icomp(nd), iparp, irtrn, j, itp
        real(c_double) xold, x, con(5*nd), rparp, y(neqp), xoutp
        logical(4) pressed
        character(1) key
        integer(c_int), parameter :: esc = 27
        common/internp/xoutp, itp

        if (nr .eq. 1) then
            itp = 1
            do j = 1, neqp
                p(j, itp) = y(j)
            end do
            xoutp = x + dz
        else
10          continue
            if (x .ge. xoutp) then
                itp = itp + 1
                do j = 1, neqp
                    p(j, itp) = contd5_p(j, xoutp, con, icomp, nd)
                end do
                xoutp = xoutp + dz
                goto 10
            end if
        end if
        return
    end subroutine soloutp

    subroutine calcpex(f, yp)
        implicit none

        real(c_double), intent(in) :: f(neqf)
        real(c_double) z, yp(neqp), artolp(1), aatolp(1), arparp(1), xoutp
        integer(c_int) i, aiparp(1), itp
        common/internp/xoutp, itp

        fp(1) = dcmplx(f(1), f(2))
        fp(2) = dcmplx(f(3), f(4))

        rparp = 0.0
        iparp = 0
        itolp = 0
        rtolp = ptol
        atolp = rtolp
        ioutp = 0
        z = zax(1)
        xoutp = z
        itp = 0
        yp = p(:, 1)
        iworkp(:) = 0
        workp(:) = 0.0d0

        artolp(1) = rtolp
        aatolp(1) = atolp
        arparp(1) = rparp
        aiparp(1) = iparp

        call dopri5_p(neqp, dpdz, z, yp, zex, artolp, aatolp, itolp, solout_fiction, 0, &
                      workp, lworkp, iworkp, liworkp, arparp, aiparp, ididp)
    end subroutine calcpex

    function dfdtcmplx(fcmplx) result(s)
        implicit none

        complex(c_double_complex) fcmplx(neqf/2), s(neqf/2), x1, x2

        x1 = xi(p(1:2*ne, :), 1)
        x2 = xi(p(2*ne + 1:4*ne, :), 1)

        s(1) = (ic*i1*x1 - ng*fcmplx(1))*q31 + 2*r1*q31*ng*cdexp(-ic*th1)*fcmplx(3) - ic*dcir1*2*q3*fcmplx(1)
        s(2) = (ic*i2*x2 - ng*fcmplx(2))*q32 + 2*r2*q32*ng*cdexp(-ic*th2)*fcmplx(3) - ic*dcir2*2*q3*fcmplx(2)
        s(3) = -fcmplx(3) + a1*fcmplx(1) + a2*fcmplx(2)

    end function dfdtcmplx

    subroutine solout_fiction
    end subroutine solout_fiction

end module fun
