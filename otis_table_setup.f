          subroutine otis_table_setup (otis_table, how_often, dt)
! Makes table of otis.f values, functions of time, with step size
! = how_often * dt

          real*8 otis_table (0:50000), dt, z, value
          integer i, j, k, how_often
          
          do i = 0, 50000
           z = dble (i) * dt * dble(how_often)
           call otis (z, value) 
           otis_table(i) = value
          end do

          end

! Time course of GABA-B, from Otis, de Koninck & Mody (1993) and proportional
! to that used in Traub et al. 1993 pyramidal cell model, J. Physiol.
                subroutine otis (t,value)

                real*8 t, value

              if (t.le.10.d0) then
                value = 0.d0
              else
            value = (1.d0 - dexp(-(t-10.d0)/38.1d0)) ** 4

       value = value * (10.2d0 * dexp(-(t-10.d0)/122.d0) +
     &    1.1d0 * dexp(-(t-10.d0)/587.d0))
              endif

                 end
