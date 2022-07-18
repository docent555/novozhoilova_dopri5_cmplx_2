program sys15f
    use, intrinsic :: iso_c_binding
    use fun
    use ifport

    implicit none

    integer(c_int) i, j, hours, minutes, seconds
    real(c_double) start_time, stop_time, calc_time
    complex(c_double_complex) pc

    call init()

    do i = 1, ne
        p(i, 1) = dreal(cdexp(ic*(i - 1)/dble(ne)*2*pi))
        p(ne + i, 1) = dimag(cdexp(ic*(i - 1)/dble(ne)*2*pi))
        p(2*ne + i, 1) = dreal(cdexp(ic*(i - 1)/dble(ne)*2*pi))
        p(3*ne + i, 1) = dimag(cdexp(ic*(i - 1)/dble(ne)*2*pi))
    end do

    write (*, '(/)')

    start_time = dclock()
    call ode4f()
    stop_time = dclock()

    calc_time = stop_time - start_time

    hours = calc_time/3600
    minutes = (calc_time - hours*3600)/60
    seconds = calc_time - hours*3600 - minutes*60

    write (*, '(/)')
    print *, 'Calcualting took:', hours, 'h :', minutes, 'm :', seconds, 's'

    do i = 1, nt
        do j = 1, 3
            fcmplx(j, i) = dcmplx(f(2*j - 1, i), f(2*j, i))
        end do
    end do

    do i = 1, nt - 1
        do j = 1, 3
            w(j, i) = imag(log(fcmplx(j, i + 1)/fcmplx(j, i)))/dt
        end do
    end do

    phi(:, 1) = 0; 
    do i = 2, nt
        do j = 1, 3
            phi(j, i) = phi(j, i - 1) + imag(log(fcmplx(j, i)/fcmplx(j, i - 1)))
        end do
    end do

    breaknum(:) = 0
    phitmp0(:) = datan2(dimag(fcmplx(:, 1)), dreal(fcmplx(:, 1)))
    phios(:, 1) = phitmp0(:)
    do i = 2, nt
        do j = 1, 3
            phitmp1(j) = datan2(dimag(fcmplx(j, i)), dreal(fcmplx(j, i)))
            if ((phitmp1(j) - phitmp0(j)) .gt. pi) breaknum(j) = breaknum(j) - 1
            if ((phitmp1(j) - phitmp0(j)) .lt. -pi) breaknum(j) = breaknum(j) + 1
            phios(j, i) = phitmp1(j) + 2.*pi*breaknum(j)
            !phios(j, i) = phitmp1(j)
            phitmp0(j) = phitmp1(j)
        end do
    end do

    do i = 1, nt - 1
        do j = 1, 3
            wos(j, i) = (phios(j, i + 1) - phios(j, i))/dt
        end do
    end do

    write (*, '(/)')

    pause

    open (1, file='F.dat')
    do i = 1, nt
        write (1, '(4e17.8)') tax(i), cdabs(fcmplx(1, i)), cdabs(fcmplx(2, i)), cdabs(fcmplx(3, i))
    end do
    close (1)

    open (1, file='FCOMP.dat')
    do i = 1, nt
        write (1, '(7e17.8)') tax(i), dreal(fcmplx(1, i)), dimag(fcmplx(1, i)), dreal(fcmplx(2, i)), &
            dimag(fcmplx(2, i)), dreal(fcmplx(3, i)), dimag(fcmplx(3, i))
    end do
    close (1)

    open (2, file='E.dat')
    do i = 1, nt
        write (2, '(5e17.8)') tax(i), eta(1, i), etag(1, i), eta(2, i), etag(2, i)
    end do
    close (2)

    open (3, file='W.dat')
    do i = 1, nt - 1
        write (3, '(4e17.8)') tax(i + 1), w(1, i), w(2, i), w(3, i)
    end do
    close (3)

    open (1, file='P.dat')
    do i = 1, nt
        write (1, '(4e17.8)') tax(i), phi(1, i), phi(2, i), phi(3, i)
    end do
    close (1)

    open (1, file='POS.dat')
    do i = 1, nt
        write (1, '(4e17.8)') tax(i), phios(1, i), phios(2, i), phios(3, i)
    end do
    close (1)

    open (3, file='WOS.dat')
    do i = 1, nt - 1
        write (3, '(4e17.8)') tax(i + 1), wos(1, i), wos(2, i), wos(3, i)
    end do
    close (3)

    stop
101 print *, 'error of file open.'
    pause
    stop
102 print *, 'error of file reading.'
    pause
    stop
103 print *, 'error of file writing.'
    pause
    stop
end program
