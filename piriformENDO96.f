c 22 March 2021: derived from piriformECT; add multipolar cells
c See /multipolar
c Multipolar cells interact with each other and L2pyr, L3pyr,
c deepbask, deepng
c For the moment, no GABA-B to multipolar
c Number of LOT input compartments for superficial interns decr.
! 13 July 2020, start piriformGJ92.f, LOT to be driven by ectopics
! 7 April 2020, start with piriformA22; for sims with gj
! 23 March 2020, from piriform156.f: try to replicate better the
! Stettler Axel 2009 result that combination of odors produce
! significant suppression.  Does this result from from stronger or
! more widespread FF inhibition?

c 23 Aug 2019, convert plateauVFOY (neocortex) to piriform cortex,
c using isomorphic structure, but different subroutines.  See notebook
c entries, 22/23 Aug 2019 and following

! 21 Dec 2018: from plateauVFOX5.f; only some tuftIB enter plateau
! 1 Dec 2018: from plateauVFO_sc8.  Generate multiple plateaus,
! with and without axonal output from tuftIB - there may be
! code around integration call, to clamp voltage of distal axons

! 23 Nov 2018: make tscale_gKDR apply to tuftIB SD only, not axon
! 19 Nov 2018: start with plateauVFO12 (20 June 2018), but make into 
! vectors of length num_tuftIB: tscale_ggabaB, tscale_gCaL &
! tscale_gKDR.  This allows plateaus to occur at different times in
! different tuftIB

! 25 May 2018, start with plateauT13, used in NIH proposal
! Then modify integrate_tuftRS to allow for VFO - see plateauAX21
! 15 15 2018, from plateau9: allow for time-varying GABA-B, gCaL, gKDR
! in tuftIB integration
! 24 Jan 2018, on Cognitive Computing Cluster
! Start with son_of_groucho133.f, use to study delta/plateau transition,
! altering GABA-B, tuftIB gCaL, etc.
! 8 March 2017, revised groucho program, start with spikewaveS96.f
! Original comments at start of that program saved in separate file ("original....")
               PROGRAM piriformA        

        PARAMETER (num_L2pyr    = 1000, 
     & num_placeholder1=100,num_placeholder2= 100,
     & num_placeholder3 = 100,
     & num_supVIP =  100, num_supng = 100, num_supintern = 500,
! supintern consists of placeholder1-3, supVIP and supng
     & num_LOT = 500, num_semilunar = 500, num_multipolar  = 500,
     & num_L3pyr = 500, 
     & num_deepbask = 200, 
     & num_deepLTS = 100, num_deepng = 200, num_deepintern = 500,
     & num_placeholder5 = 500,
     & num_placeholder6 = 500)

        PARAMETER (ncellspernode = 500, 
     &    nodesfor_L2pyr     = 2,
     &    nodesfor_supintern   = 1,
     &    nodesfor_LOT       = 1,
     &    nodesfor_semilunar = 1,
     &    nodesfor_multipolar   = 1,
     &    nodesfor_L3pyr     = 1,
     &    nodesfor_deepintern  = 1,
     &    nodesfor_placeholder5 = 1,
     &    nodesfor_placeholder6 = 1)

        PARAMETER (numnodes = 10)  ! Check manually for consistency.
        PARAMETER (maxcellspernode = 500)
            

        PARAMETER (numcomp_L2pyr      = 74,
     &             numcomp_supVIP     = 59,
     &             numcomp_placeholder1 = 59,
     &             numcomp_placeholder2 = 59,
     &             numcomp_placeholder3 = 59,
     &             numcomp_LOT          = 59,
     &             numcomp_semilunar    = 74,
     &             numcomp_multipolar   = 59,
     &             numcomp_L3pyr        = 74,
     &             numcomp_deepbask     = 59,
     &             numcomp_deepLTS      = 59,
     &             numcomp_placeholder5 =137,
     &             numcomp_placeholder6 = 74,
     &             numcomp_supng        = 59,
     &             numcomp_deepng       = 59)

c       PARAMETER (num_L2pyr_to_L2pyr = 50,
c       PARAMETER (num_L2pyr_to_L2pyr = 15,
        PARAMETER (num_L2pyr_to_L2pyr = 20,
c       PARAMETER (num_L2pyr_to_L2pyr =  5,
!       PARAMETER (num_L2pyr_to_L2pyr = 10,
     &   num_L2pyr_to_placeholder1   =  1,
     &   num_L2pyr_to_placeholder2   =  1,
     &   num_L2pyr_to_placeholder3    = 1, 
     &   num_L2pyr_to_supng     = 10, ! note
     &   num_L2pyr_to_LOT =  1, ! make sure conductance = 0 
     &   num_L2pyr_to_semilunar    =  2,
c    &   num_L2pyr_to_multipolar      = 20,
     &   num_L2pyr_to_multipolar      =  5,
     &   num_L2pyr_to_deepbask  = 60,
     &   num_L2pyr_to_deepLTS  = 60,
     &   num_L2pyr_to_deepng   = 60,
     &   num_L2pyr_to_supVIP    = 10, 
! No L2pyr to deepng - in original program, change this
c    &   num_L2pyr_to_L3pyr = 50)
c    &   num_L2pyr_to_L3pyr =  5)
     &   num_L2pyr_to_L3pyr = 10)

        PARAMETER
     &  (num_placeholder1_to_L2pyr   =  1,
     &   num_placeholder1_to_placeholder1    =  1,
     &   num_placeholder1_to_placeholder2    =  1,
     &   num_placeholder1_to_placeholder3     =  1,
     &   num_placeholder1_to_supng      =  1,
     &   num_placeholder1_to_LOT  =  1,

     &   num_placeholder2_to_L2pyr   =  1, ! note
     &   num_placeholder2_to_LOT  = 1,
     &   num_placeholder2_to_semilunar     = 1,
     &   num_placeholder2_to_multipolar       = 1,
     &   num_placeholder2_to_L3pyr  = 5,

     &   num_placeholder3_to_placeholder1     =  1,
     &   num_placeholder3_to_placeholder2     =  1,
     &   num_placeholder3_to_placeholder3      =  1,
     &   num_placeholder3_to_L2pyr    =  1,
     &   num_placeholder3_to_LOT   =  1,
     &   num_placeholder3_to_semilunar      =  1)

        PARAMETER
     &  (num_supng_to_L2pyr     = 20,
     &   num_supng_to_L3pyr    = 20,
     &   num_supng_to_semilunar       = 20,
     &   num_supng_to_multipolar         =  1,
     &   num_supng_to_supng        =  5,
     &   num_supng_to_placeholder1      = 1)

        PARAMETER
     &  (num_placeholder3_to_multipolar        =  1,
     &   num_placeholder3_to_deepbask    =  1,
     &   num_placeholder3_to_deepLTS    =  1,
     &   num_placeholder3_to_supVIP      =  1,
     &   num_placeholder3_to_L3pyr   =  1,

     &   num_LOT_to_L2pyr = 20,
     &   num_LOT_to_placeholder1  =  1,
     &   num_LOT_to_placeholder2  =  1,
     &   num_LOT_to_placeholder3   =  1,
     &   num_LOT_to_LOT=  1,
     &   num_LOT_to_semilunar   = 30,
     &   num_LOT_to_multipolar     =  1,
     &   num_LOT_to_deepbask =  1,
     &   num_LOT_to_deepLTS =  1,
     &   num_LOT_to_supVIP   = 30,
     &   num_LOT_to_supng    = 30,
     &   num_LOT_to_deepng   =  1,
     &   num_LOT_to_L3pyr= 30,

     &   num_semilunar_to_L2pyr    = 40,
     &   num_semilunar_to_placeholder1     = 1) 

        PARAMETER
     &  (num_semilunar_to_placeholder2     = 1,
     &   num_semilunar_to_placeholder3      = 1, 
     &   num_semilunar_to_LOT   =  1,
     &   num_semilunar_to_semilunar      = 2,
     &   num_semilunar_to_multipolar        =  1,
     &   num_semilunar_to_deepbask    = 20,
     &   num_semilunar_to_deepLTS    = 20,
     &   num_semilunar_to_supVIP      = 5,
     &   num_semilunar_to_deepng      = 20,
     &   num_semilunar_to_L3pyr   = 40,
c ? include semilunar to supng?

     &   num_multipolar_to_L2pyr    = 15,
     &   num_multipolar_to_placeholder1     =  1,
     &   num_multipolar_to_placeholder2     =  1,
     &   num_multipolar_to_placeholder3      =  1,
     &   num_multipolar_to_LOT   =  1,
     &   num_multipolar_to_semilunar      =  1,
c    &   num_multipolar_to_multipolar      = 25,
     &   num_multipolar_to_multipolar      =  3,
     &   num_multipolar_to_deepbask    = 10,
     &   num_multipolar_to_deepLTS    =  1,
     &   num_multipolar_to_supVIP      =  1,
     &   num_multipolar_to_deepng      = 10,
     &   num_multipolar_to_L3pyr   = 15)

        PARAMETER
     &  (num_deepbask_to_LOT =  1,
     &   num_deepbask_to_semilunar    = 20,
     &   num_deepbask_to_multipolar      = 15,
     &   num_deepbask_to_deepbask  = 10,
     &   num_deepbask_to_deepLTS  = 10,
     &   num_deepbask_to_supVIP    =  1,
     &   num_deepbask_to_deepng    = 10,
     &   num_deepbask_to_L2pyr = 20,
     &   num_deepbask_to_L3pyr = 20,

     &   num_deepLTS_to_L2pyr  = 20,
     &   num_deepLTS_to_LOT =  1,
     &   num_deepLTS_to_semilunar    = 20,
     &   num_deepLTS_to_multipolar      =  1,
     &   num_deepLTS_to_L3pyr = 20,

     &   num_supVIP_to_L2pyr   = 20)
        PARAMETER
     &  (
     &   num_supVIP_to_placeholder1    =  1,
     &   num_supVIP_to_placeholder2    =  1,
     &   num_supVIP_to_placeholder3     =  1,
     &   num_supVIP_to_LOT  =  1,
     &   num_supVIP_to_semilunar     = 20,
     &   num_supVIP_to_multipolar       =  1,
     &   num_supVIP_to_deepbask   =  1,
     &   num_supVIP_to_deepLTS   =  1,
     &   num_supVIP_to_supVIP    =  1,
     &   num_supVIP_to_supng     =   1,
     &   num_supVIP_to_L3pyr  = 20)

        PARAMETER
     &  (num_deepng_to_semilunar      = 20,
     &   num_deepng_to_multipolar        =  1,
     &   num_deepng_to_L2pyr   = 20,
     &   num_deepng_to_L3pyr   = 20,
     &   num_deepng_to_LOT   =  1,
     &   num_deepng_to_deepng      =  2,
     &   num_deepng_to_deepbask    =  2,

     &   num_placeholder5_to_L2pyr       =  1,
     &   num_placeholder5_to_placeholder1        =  1,
     &   num_placeholder5_to_placeholder2        =  1,
     &   num_placeholder5_to_supng          =  1,
     &   num_placeholder5_to_LOT      =  1,
     &   num_placeholder5_to_semilunar         =  1,
     &   num_placeholder5_to_multipolar           =  1,
     &   num_placeholder5_to_deepbask       =  1,
     &   num_placeholder5_to_deepLTS       =  1,
     &   num_placeholder5_to_deepng         =  1,
     &   num_placeholder5_to_placeholder6            =  1,       
     &   num_placeholder5_to_L3pyr      =  1,

     &   num_placeholder6_to_placeholder5            =  1, ! note
     &   num_placeholder6_to_placeholder6            =  1)
        PARAMETER
c    &  (num_L3pyr_to_L2pyr = 20,
c    &  (num_L3pyr_to_L2pyr =  5,
     &  (num_L3pyr_to_L2pyr = 10,
     &   num_L3pyr_to_placeholder1  =  1,
     &   num_L3pyr_to_placeholder2  =  1,
     &   num_L3pyr_to_placeholder3   =  1,
     &   num_L3pyr_to_LOT=  1,
     &   num_L3pyr_to_semilunar   =  2,
c    &   num_L3pyr_to_multipolar     = 20,
     &   num_L3pyr_to_multipolar     =  5,
     &   num_L3pyr_to_deepbask = 40,
     &   num_L3pyr_to_deepLTS = 40,
     &   num_L3pyr_to_supVIP   =  1,
     &   num_L3pyr_to_deepng   = 40,
     &   num_L3pyr_to_placeholder5      =  1,
     &   num_L3pyr_to_placeholder6      =  1,
     &   num_L3pyr_to_L3pyr= 10)
c    &   num_L3pyr_to_L3pyr=  5)
c    &   num_L3pyr_to_L3pyr= 20)

c Begin definition of number of compartments that can be
c contacted for each type of synaptic connection.
        PARAMETER (ncompallow_L2pyr_to_L2pyr = 36,
     &   ncompallow_L2pyr_to_placeholder1   =  1,
     &   ncompallow_L2pyr_to_placeholder2   =  1,
     &   ncompallow_L2pyr_to_placeholder3    =  1,
     &   ncompallow_L2pyr_to_supng     = 52,
     &   ncompallow_L2pyr_to_LOT =  1,
     &   ncompallow_L2pyr_to_semilunar    =  7,
     &   ncompallow_L2pyr_to_multipolar      = 24,
     &   ncompallow_L2pyr_to_deepbask  = 24,
     &   ncompallow_L2pyr_to_deepLTS  = 24,
     &   ncompallow_L2pyr_to_deepng   = 52,
     &   ncompallow_L2pyr_to_supVIP    = 24,
     &   ncompallow_L2pyr_to_L3pyr = 36)

        PARAMETER (ncompallow_placeholder1_to_L2pyr   =  1,
     &   ncompallow_placeholder1_to_placeholder1     =  1,
     &   ncompallow_placeholder1_to_supng       = 1, 
     &   ncompallow_placeholder1_to_placeholder2     =  1,
     &   ncompallow_placeholder1_to_placeholder3      =  1,
     &   ncompallow_placeholder1_to_LOT   =  1)

        PARAMETER (ncompallow_placeholder3_to_L2pyr    =  1,
     &   ncompallow_placeholder3_to_placeholder1      =  1,
     &   ncompallow_placeholder3_to_placeholder2      =  1,
     &   ncompallow_placeholder3_to_placeholder3       =  1,
     &   ncompallow_placeholder3_to_LOT    =  1,
     &   ncompallow_placeholder3_to_semilunar       =  1,
     &   ncompallow_placeholder3_to_multipolar         =  1,
     &   ncompallow_placeholder3_to_deepbask     =  1,
     &   ncompallow_placeholder3_to_deepLTS     =  1,
     &   ncompallow_placeholder3_to_supVIP       =  1,
     &   ncompallow_placeholder3_to_L3pyr    =  1)

        PARAMETER (ncompallow_supng_to_L2pyr = 24,
     &   ncompallow_supng_to_L3pyr     = 24,
     &   ncompallow_supng_to_semilunar        = 24,
     &   ncompallow_supng_to_multipolar          =  1,
     &   ncompallow_supng_to_supng         =  4,
     &   ncompallow_supng_to_placeholder1       =  1)

        PARAMETER (ncompallow_LOT_to_L2pyr = 24,
     &   ncompallow_LOT_to_placeholder1   =  1,
     &   ncompallow_LOT_to_placeholder2   =  1,
     &   ncompallow_LOT_to_placeholder3    =  1,
     &   ncompallow_LOT_to_LOT =  1,
     &   ncompallow_LOT_to_semilunar    = 24,
     &   ncompallow_LOT_to_multipolar      =  1,
     &   ncompallow_LOT_to_deepbask  =  1,
     &   ncompallow_LOT_to_deepLTS  =  1,
c    &   ncompallow_LOT_to_supVIP    = 24,
     &   ncompallow_LOT_to_supVIP    =  9,
c    &   ncompallow_LOT_to_supng     = 24,
     &   ncompallow_LOT_to_supng     =  9, ! corrected 22Mar2021
     &   ncompallow_LOT_to_deepng    =  1,
     &   ncompallow_LOT_to_L3pyr = 24)

        PARAMETER (ncompallow_semilunar_to_L2pyr   = 43,
     &   ncompallow_semilunar_to_placeholder1      =  1,
     &   ncompallow_semilunar_to_placeholder2      =  1,
     &   ncompallow_semilunar_to_placeholder3       =  1,
     &   ncompallow_semilunar_to_LOT    =  1,
     &   ncompallow_semilunar_to_semilunar       = 15,
     &   ncompallow_semilunar_to_multipolar         =  1,
     &   ncompallow_semilunar_to_deepbask     = 24,
     &   ncompallow_semilunar_to_deepLTS     = 24,
     &   ncompallow_semilunar_to_supVIP       = 24,
     &   ncompallow_semilunar_to_deepng       = 52,
     &   ncompallow_semilunar_to_L3pyr    = 43)

        PARAMETER (ncompallow_multipolar_to_L2pyr   = 36,
     &   ncompallow_multipolar_to_placeholder1      =  1,
     &   ncompallow_multipolar_to_placeholder2      =  1,
     &   ncompallow_multipolar_to_placeholder3       =  1,
     &   ncompallow_multipolar_to_LOT    =  1,
     &   ncompallow_multipolar_to_semilunar       =  1,
     &   ncompallow_multipolar_to_multipolar       = 24,
     &   ncompallow_multipolar_to_deepbask     = 24,
     &   ncompallow_multipolar_to_deepLTS     =  1,
     &   ncompallow_multipolar_to_supVIP       =  1,
     &   ncompallow_multipolar_to_deepng       = 24,
     &   ncompallow_multipolar_to_L3pyr    = 36)

        PARAMETER (ncompallow_deepbask_to_LOT = 1, 
     &   ncompallow_deepbask_to_semilunar     = 11,
     &   ncompallow_deepbask_to_multipolar       = 24,
     &   ncompallow_deepbask_to_deepbask   = 24,
     &   ncompallow_deepbask_to_deepLTS   = 24,
     &   ncompallow_deepbask_to_supVIP     = 24,
     &   ncompallow_deepbask_to_deepng     =  4,
     &   ncompallow_deepbask_to_L2pyr  = 11,
     &   ncompallow_deepbask_to_L3pyr  = 11)

        PARAMETER (ncompallow_supVIP_to_L2pyr = 24,
     &   ncompallow_supVIP_to_placeholder1     =  1,
     &   ncompallow_supVIP_to_placeholder2     =  1,
     &   ncompallow_supVIP_to_placeholder3      =  1,
     &   ncompallow_supVIP_to_supng       = 20,
     &   ncompallow_supVIP_to_LOT   =  1,
     &   ncompallow_supVIP_to_semilunar      = 24, ! should equal the number of VIP cells to each semilunar
     &   ncompallow_supVIP_to_multipolar        =  1,
     &   ncompallow_supVIP_to_deepbask    =  4,
     &   ncompallow_supVIP_to_deepLTS    =  4,
     &   ncompallow_supVIP_to_supVIP     =  4,
     &   ncompallow_supVIP_to_L3pyr   = 24)

        PARAMETER (ncompallow_deepng_to_semilunar = 33,
     &   ncompallow_deepng_to_multipolar      = 24,
     &   ncompallow_deepng_to_L2pyr = 33,
     &   ncompallow_deepng_to_L3pyr = 33,
     &   ncompallow_deepng_to_LOT =  1,
     &   ncompallow_deepng_to_deepng    =  4,
     &   ncompallow_deepng_to_deepbask  = 52)

        PARAMETER (ncompallow_deepLTS_to_L2pyr = 24,
     &   ncompallow_deepLTS_to_LOT = 1,
     &   ncompallow_deepLTS_to_semilunar = 24,
     &   ncompallow_deepLTS_to_multipolar   = 1,
     &   ncompallow_deepLTS_to_L3pyr = 24)

        PARAMETER (ncompallow_placeholder5_to_L2pyr =  1,
     &   ncompallow_placeholder5_to_placeholder1      =  1,
     &   ncompallow_placeholder5_to_placeholder2      =  1,
     &   ncompallow_placeholder5_to_supng        =  1,
     &   ncompallow_placeholder5_to_LOT    =  1,
     &   ncompallow_placeholder5_to_semilunar       =  1,
     &   ncompallow_placeholder5_to_multipolar         =  1,
     &   ncompallow_placeholder5_to_deepbask     =  1,
!    &   ncompallow_placeholder5_to_deepLTS     =  1,
     &   ncompallow_placeholder5_to_deepLTS     =  1,
     &   ncompallow_placeholder5_to_deepng       =  1, 
     &   ncompallow_placeholder5_to_placeholder6          =  1,
     &   ncompallow_placeholder5_to_L3pyr    =  1)

        PARAMETER (ncompallow_placeholder6_to_placeholder5 =  1,
     &   ncompallow_placeholder6_to_placeholder6 =  1)

        PARAMETER (ncompallow_L3pyr_to_L2pyr = 36,
     &    ncompallow_L3pyr_to_placeholder1   =  1,
     &    ncompallow_L3pyr_to_placeholder2   =  1,
     &    ncompallow_L3pyr_to_placeholder3    =  1,
     &    ncompallow_L3pyr_to_LOT =  1,
     &    ncompallow_L3pyr_to_semilunar    =  7,
     &    ncompallow_L3pyr_to_multipolar      = 24,
     &    ncompallow_L3pyr_to_deepbask  = 24,
     &    ncompallow_L3pyr_to_deepLTS  = 24,
     &    ncompallow_L3pyr_to_supVIP    = 24,
     &    ncompallow_L3pyr_to_deepng    = 52,
     &    ncompallow_L3pyr_to_placeholder5       = 1,
     &    ncompallow_L3pyr_to_placeholder6       =  1,
     &    ncompallow_L3pyr_to_L3pyr = 36)
c End definition of number of allowed compartments that
c can be contacted for each sort of connection

c Note that gj form only between cells of a given type and in same node
c Except different sorts of interneurons may couple
c gj/cell = 2 x total gj / # cells
c for proportions, see /home/traub/supergj/tests.f
c???? GJ BETWEEN SUPVIP ????
       integer, parameter :: totaxgj_L2pyr =  100
c      integer, parameter :: totaxgj_L2pyr =    1
c      integer, parameter :: totaxgj_L2pyr = 1000
c      integer, parameter :: totaxgj_L2pyr = 1200
       integer, parameter :: totSDgj_placeholder1   =   1 
       integer, parameter :: totSDgj_placeholder2   =   1 
       integer, parameter :: totSDgj_placeholder3    =  1  
       integer, parameter :: totaxgj_LOT =   1 
c      integer, parameter :: totaxgj_semilunar    = 750 
       integer, parameter :: totaxgj_semilunar    = 200 
       integer, parameter :: totaxgj_multipolar   = 500 ! does code
c allow for axonal gj?
       integer, parameter :: totaxgj_L3pyr = 200 
       integer, parameter :: totSDgj_deepbask  = 250 
       integer, parameter :: totSDgj_deepLTS  =   1 
c      integer, parameter :: totSDgj_supVIP    = 250 
       integer, parameter :: totSDgj_supVIP    =1000 
       integer, parameter :: totaxgj_placeholder5       =   1  
       integer, parameter :: totaxgj_placeholder6       =   1
       integer, parameter :: totSDgj_supng     = 250
       integer, parameter :: totSDgj_deepng    = 250

c Define number of compartments on a cell where a gj might form
       integer, parameter :: num_axgjcompallow_L2pyr = 1
       integer, parameter :: num_SDgjcompallow_placeholder1  = 8
       integer, parameter :: num_SDgjcompallow_supng    = 8
       integer, parameter :: num_SDgjcompallow_placeholder3   = 8
       integer, parameter :: num_axgjcompallow_LOT= 1
       integer, parameter :: num_axgjcompallow_semilunar   = 1
       integer, parameter :: num_axgjcompallow_multipolar     = 1
       integer, parameter :: num_axgjcompallow_L3pyr= 1
       integer, parameter :: num_SDgjcompallow_deepbask = 8
       integer, parameter :: num_SDgjcompallow_deepng   = 8
       integer, parameter :: num_SDgjcompallow_supVIP   = 8
       integer, parameter :: num_axgjcompallow_placeholder5    = 1
       integer, parameter :: num_axgjcompallow_placeholder6    = 1

c Define gap junction conductances.
!      double precision, parameter :: gapcon_L2pyr  = 6.0d-3 ! also
!   define as just double precision, so as to be able to vary it
       double precision, parameter :: gapcon_placeholder1   = 0.d-3
       double precision, parameter :: gapcon_supng     = 0.5d-3
       double precision, parameter :: gapcon_placeholder2   = 0.d-3
       double precision, parameter :: gapcon_placeholder3    = 0.d-3

       double precision, parameter :: gapcon_LOT = 0.d-3 

       double precision, parameter :: gapcon_semilunar    = 3.00d-3
       double precision, parameter :: gapcon_multipolar      = 0.d-3
       double precision, parameter :: gapcon_L3pyr = 3.00d-3
c      double precision, parameter :: gapcon_L3pyr = 0.d-3
c      double precision, parameter :: gapcon_L3pyr = 8.d-3
       double precision, parameter :: gapcon_deepbask  = 1.d-3
       double precision, parameter :: gapcon_deepng    = 0.5d-3
       double precision, parameter :: gapcon_deepLTS  = 0.d-3
c      double precision, parameter :: gapcon_supVIP    = 2.d-3
       double precision, parameter :: gapcon_supVIP    = 0.1d-3
       double precision, parameter :: gapcon_placeholder5    = 0.d-3
       double precision, parameter :: gapcon_placeholder6    = 0.0d-3


c Assorted parameters
         double precision, parameter :: dt = 0.002d0
         double precision, parameter :: Mg = 1.00 ! for NMDA-dependent CCh delta, try lower Mg
! Castro-Alamancos J Physiol, disinhib. neocortex in vitro, uses
! Mg = 1.3
         double precision, parameter :: NMDA_saturation_fact
!    &                                   = 5.d0
     &                                   = 80.d0
c NMDA conductance developed on one postsynaptic compartment,
c from one type of presynaptic cell, can be at most this
c factor x unitary conductance
c UNFORTUNATELY, with this scheme,if one NMDA cond. set to 0
c on a cell type, all NMDA conductances will be forced to 0
c on that cell type...

       double precision, parameter :: thal_cort_delay = 1.d0
       double precision, parameter :: cort_thal_delay = 5.d0
       integer, parameter :: how_often = 50
! how_often defines how many time steps between synaptic conductance
! updates, and between broadcastings of axonal voltages.
       double precision, parameter :: axon_refrac_time = 1.5d0

c For these ectopic rate parameters, assume noisepe checked
c every 200 time steps = 0.4 ms = 1./2.5 ms
      double precision, parameter :: noisepe_L2pyr   =
     &      0.d0 
c    &            1.d0 / (2.5d0 * 1000.d0)
c    &            1.d0 / (2.5d0 * 2000.d0)
      double precision, parameter :: noisepe_LOT  =
c    &            1.d0 / (2.5d0 *  300.d0)
     &            1.d0 / (2.5d0 *  400.d0)
c    &            0.d0 
! Note that noisepe_semilunar will be time-dependent
      double precision, parameter :: noisepe_semilunar_save     =
     &            0.d0 / (2.5d0 * 5000.d0)
      double precision, parameter :: noisepe_multipolar_save =
c this one also will be time_dependent
     &            0.d0 / (2.5d0 * 800.d0)
      double precision, parameter :: noisepe_L3pyr  =
c    &            1.d0 / (2.5d0 * 1000.d0)
     &            0.d0 / (2.5d0 * 1000.d0)
      double precision, parameter :: noisepe_placeholder5        =
     &            0.d0 / (2.5d0 * 20000.d0)


c Synaptic conductance time constants. 
      real*8, parameter :: tauAMPA_L2pyr_to_L2pyr=2.d0 
      real*8, parameter :: tauNMDA_L2pyr_to_L2pyr=130.5d0 
      real*8, parameter :: tauAMPA_L2pyr_to_placeholder1  =.8d0   
      real*8, parameter :: tauNMDA_L2pyr_to_placeholder1  =100.d0 
      real*8, parameter :: tauAMPA_L2pyr_to_supng    =.8d0   
      real*8, parameter :: tauNMDA_L2pyr_to_supng    =100.d0 
      real*8, parameter :: tauAMPA_L2pyr_to_placeholder2  =.8d0  
      real*8, parameter :: tauNMDA_L2pyr_to_placeholder2  =100.d0 
      real*8, parameter :: tauAMPA_L2pyr_to_placeholder3   =1.d0  
      real*8, parameter :: tauNMDA_L2pyr_to_placeholder3   =100.d0 
      real*8, parameter :: tauAMPA_L2pyr_to_LOT=2.d0   
      real*8, parameter :: tauNMDA_L2pyr_to_LOT=130.d0 
      real*8, parameter :: tauAMPA_L2pyr_to_semilunar   =2.d0   
      real*8, parameter :: tauNMDA_L2pyr_to_semilunar   =130.d0 
      real*8, parameter :: tauAMPA_L2pyr_to_multipolar     =2.d0   
      real*8, parameter :: tauNMDA_L2pyr_to_multipolar   =130.d0 
      real*8, parameter :: tauAMPA_L2pyr_to_deepbask =.8d0   
      real*8, parameter :: tauNMDA_L2pyr_to_deepbask =100.d0 
      real*8, parameter :: tauAMPA_L2pyr_to_deepng   =.8d0   
      real*8, parameter :: tauNMDA_L2pyr_to_deepng   =100.d0 
      real*8, parameter :: tauAMPA_L2pyr_to_deepLTS =.8d0   
      real*8, parameter :: tauNMDA_L2pyr_to_deepLTS =100.d0 
      real*8, parameter :: tauAMPA_L2pyr_to_supVIP   =1.d0   
      real*8, parameter :: tauNMDA_L2pyr_to_supVIP   =100.d0 
      real*8, parameter :: tauAMPA_L2pyr_to_L3pyr=2.d0   
      real*8, parameter :: tauNMDA_L2pyr_to_L3pyr=130.d0 

      real*8,  parameter :: tauGABA_placeholder1_to_L2pyr   =6.d0
      real*8,  parameter :: tauGABA_placeholder1_to_placeholder1 =3.d0  
      real*8,  parameter :: tauGABA_placeholder1_to_supng      =3.d0  
      real*8,  parameter :: tauGABA_placeholder1_to_placeholder2 =3.d0  
      real*8,  parameter :: tauGABA_placeholder1_to_placeholder3 =3.d0  
      real*8,  parameter :: tauGABA_placeholder1_to_LOT  =6.d0 

      real*8,  parameter :: tauGABA_placeholder2_to_L2pyr   =6.d0 
      real*8,  parameter :: tauGABA_placeholder2_to_LOT  =6.d0 
      real*8,  parameter :: tauGABA_placeholder2_to_semilunar    =6.d0 
      real*8,  parameter :: tauGABA_placeholder2_to_multipolar =6.d0 
      real*8,  parameter :: tauGABA_placeholder2_to_L3pyr  =6.d0 

      real*8, parameter :: tauGABA_placeholder3_to_L2pyr    =20.d0 
      real*8, parameter :: tauGABA_placeholder3_to_placeholder1 =20.d0 
      real*8, parameter :: tauGABA_placeholder3_to_placeholder2 =20.d0 
      real*8, parameter :: tauGABA_placeholder3_to_placeholder3 =20.d0 
      real*8, parameter :: tauGABA_placeholder3_to_LOT   =20.d0 
      real*8, parameter :: tauGABA_placeholder3_to_semilunar     =20.d0 
      real*8, parameter :: tauGABA_placeholder3_to_multipolar  =20.d0 
      real*8, parameter :: tauGABA_placeholder3_to_deepbask    =20.d0 
      real*8, parameter :: tauGABA_placeholder3_to_deepLTS    =20.d0 
      real*8, parameter :: tauGABA_placeholder3_to_supVIP      =20.d0 
      real*8, parameter :: tauGABA_placeholder3_to_L3pyr   =20.d0  

      real*8, parameter:: tauGABA_supng_to_L2pyr      =6.d0
      real*8, parameter:: tauGABAB_supng_to_L2pyr     =100.d0
      real*8, parameter:: tauGABA_supng_to_L3pyr     =6.d0
      real*8, parameter:: tauGABAB_supng_to_L3pyr    =100.d0
      real*8, parameter:: tauGABA_supng_to_semilunar        =6.d0
      real*8, parameter:: tauGABAB_supng_to_semilunar       =100.d0
      real*8, parameter:: tauGABA_supng_to_multipolar        =6.d0
      real*8, parameter:: tauGABAB_supng_to_multipolar       =100.d0
      real*8, parameter:: tauGABA_supng_to_supng         =3.d0
      real*8, parameter:: tauGABA_supng_to_placeholder1       =3.d0

      real*8, parameter :: tauAMPA_LOT_to_L2pyr =2.d0  
      real*8, parameter :: tauNMDA_LOT_to_L2pyr =130.d0 
      real*8, parameter :: tauAMPA_LOT_to_placeholder1  =.8d0  
      real*8, parameter :: tauNMDA_LOT_to_placeholder1  =100.d0
      real*8, parameter :: tauAMPA_LOT_to_placeholder2  =.8d0  
      real*8, parameter :: tauNMDA_LOT_to_placeholder2  =100.d0
      real*8, parameter :: tauAMPA_LOT_to_placeholder3   =1.d0  
      real*8, parameter :: tauNMDA_LOT_to_placeholder3   =100.d0
      real*8, parameter :: tauAMPA_LOT_to_LOT=2.d0  
!     real*8, parameter :: tauNMDA_LOT_to_LOT=130.d0 
      real*8, parameter :: tauNMDA_LOT_to_LOT= 15.d0 ! small tau per Fleidervish et al., NEURON 
      real*8, parameter :: tauAMPA_LOT_to_semilunar   =2.d0  
      real*8, parameter :: tauNMDA_LOT_to_semilunar   =130.d0 
      real*8, parameter :: tauAMPA_LOT_to_multipolar   =2.d0  
      real*8, parameter :: tauNMDA_LOT_to_multipolar   =130.d0
      real*8, parameter :: tauAMPA_LOT_to_deepbask =.8d0  
      real*8, parameter :: tauNMDA_LOT_to_deepbask =100.d0
      real*8, parameter :: tauAMPA_LOT_to_deepng   =.8d0  
      real*8, parameter :: tauNMDA_LOT_to_deepng   =100.d0
      real*8, parameter :: tauAMPA_LOT_to_deepLTS =.8d0  
      real*8, parameter :: tauNMDA_LOT_to_deepLTS =100.d0
      real*8, parameter :: tauAMPA_LOT_to_supVIP   =1.d0  
      real*8, parameter :: tauNMDA_LOT_to_supVIP   =100.d0
      real*8, parameter :: tauAMPA_LOT_to_supng    =1.d0  
      real*8, parameter :: tauNMDA_LOT_to_supng    =100.d0
      real*8, parameter :: tauAMPA_LOT_to_L3pyr=2.d0  
      real*8, parameter :: tauNMDA_LOT_to_L3pyr=130.d0

      real*8, parameter :: tauAMPA_semilunar_to_L2pyr    =2.d0 
      real*8, parameter :: tauNMDA_semilunar_to_L2pyr    =130.d0
      real*8, parameter :: tauAMPA_semilunar_to_placeholder1     =.8d0  
      real*8, parameter :: tauNMDA_semilunar_to_placeholder1    =100.d0 
      real*8, parameter :: tauAMPA_semilunar_to_placeholder2     =.8d0  
      real*8, parameter :: tauNMDA_semilunar_to_placeholder2    =100.d0 
      real*8, parameter :: tauAMPA_semilunar_to_placeholder3      =1.d0  
      real*8, parameter :: tauNMDA_semilunar_to_placeholder3    =100.d0 
      real*8, parameter :: tauAMPA_semilunar_to_LOT   =2.d0   
      real*8, parameter :: tauNMDA_semilunar_to_LOT   =130.d0 
      real*8, parameter :: tauAMPA_semilunar_to_semilunar      =2.d0  
      real*8, parameter :: tauNMDA_semilunar_to_semilunar      =130.d0 
      real*8, parameter :: tauAMPA_semilunar_to_multipolar     =2.0d0 
      real*8, parameter :: tauNMDA_semilunar_to_multipolar    =130.d0 
      real*8, parameter :: tauAMPA_semilunar_to_deepbask    =.8d0  
      real*8, parameter :: tauNMDA_semilunar_to_deepbask    =100.d0 
      real*8, parameter :: tauAMPA_semilunar_to_deepng      =.8d0  
      real*8, parameter :: tauNMDA_semilunar_to_deepng      =100.d0 
      real*8, parameter :: tauAMPA_semilunar_to_deepLTS    =.8d0  
      real*8, parameter :: tauNMDA_semilunar_to_deepLTS    =100.d0 
      real*8, parameter :: tauAMPA_semilunar_to_supVIP      =1.d0  
      real*8, parameter :: tauNMDA_semilunar_to_supVIP      =100.d0 
      real*8, parameter :: tauAMPA_semilunar_to_L3pyr   =2.0d0 
      real*8, parameter :: tauNMDA_semilunar_to_L3pyr   =130.d0 

      real*8, parameter :: tauAMPA_multipolar_to_L2pyr    =2.d0 
      real*8, parameter :: tauNMDA_multipolar_to_L2pyr    =130.d0
      real*8, parameter :: tauAMPA_multipolar_to_placeholder1   =.8d0  
      real*8, parameter :: tauNMDA_multipolar_to_placeholder1 =100.d0 
      real*8, parameter :: tauAMPA_multipolar_to_placeholder2   =.8d0  
      real*8, parameter :: tauNMDA_multipolar_to_placeholder2 =100.d0 
      real*8, parameter :: tauAMPA_multipolar_to_placeholder3   =1.d0  
      real*8, parameter :: tauNMDA_multipolar_to_placeholder3 =100.d0 
      real*8, parameter :: tauAMPA_multipolar_to_LOT   =2.d0  
      real*8, parameter :: tauNMDA_multipolar_to_LOT   =130.d0 
      real*8, parameter :: tauAMPA_multipolar_to_semilunar     =2.d0  
      real*8, parameter :: tauNMDA_multipolar_to_semilunar    =130.d0 
      real*8, parameter :: tauAMPA_multipolar_to_multipolar     =2.d0  
      real*8, parameter :: tauNMDA_multipolar_to_multipolar   =130.d0 
      real*8, parameter :: tauAMPA_multipolar_to_deepbask    =.8d0  
      real*8, parameter :: tauNMDA_multipolar_to_deepbask    =100.d0 
      real*8, parameter :: tauAMPA_multipolar_to_deepng      =.8d0  
      real*8, parameter :: tauNMDA_multipolar_to_deepng      =100.d0 
      real*8, parameter :: tauAMPA_multipolar_to_deepLTS    =.8d0  
      real*8, parameter :: tauNMDA_multipolar_to_deepLTS    =100.d0 
      real*8, parameter :: tauAMPA_multipolar_to_supVIP      =1.d0   
      real*8, parameter :: tauNMDA_multipolar_to_supVIP      =100.d0 
      real*8, parameter :: tauAMPA_multipolar_to_L3pyr   =2.d0  
      real*8, parameter :: tauNMDA_multipolar_to_L3pyr   =130.d0 

      real*8, parameter :: tauGABA_deepbask_to_LOT =6.d0 
      real*8, parameter :: tauGABA_deepbask_to_semilunar    =6.d0  
      real*8, parameter :: tauGABA_deepbask_to_multipolar      =6.d0  
      real*8, parameter :: tauGABA_deepbask_to_deepbask  =3.d0  
      real*8, parameter :: tauGABA_deepbask_to_deepLTS  =3.d0  
      real*8, parameter :: tauGABA_deepbask_to_supVIP    =3.d0  
      real*8, parameter :: tauGABA_deepbask_to_deepng    =3.d0  
      real*8, parameter :: tauGABA_deepbask_to_L2pyr =6.d0  
      real*8, parameter :: tauGABA_deepbask_to_L3pyr =6.d0  

      real*8, parameter :: tauGABA_deepLTS_to_L2pyr   =6.d0 
      real*8, parameter :: tauGABA_deepLTS_to_LOT  =6.d0 
      real*8, parameter :: tauGABA_deepLTS_to_semilunar   =6.d0 
      real*8, parameter :: tauGABA_deepLTS_to_multipolar=6.d0 
      real*8, parameter :: tauGABA_deepLTS_to_L3pyr  =6.d0 

      real*8, parameter :: tauGABA_supVIP_to_L2pyr    =20.d0 
      real*8, parameter :: tauGABA_supVIP_to_placeholder1 =20.d0 
      real*8, parameter :: tauGABA_supVIP_to_placeholder2 =20.d0 
      real*8, parameter :: tauGABA_supVIP_to_placeholder3 =20.d0 
      real*8, parameter :: tauGABA_supVIP_to_supng       =20.d0 
      real*8, parameter :: tauGABA_supVIP_to_LOT   =20.d0 
      real*8, parameter :: tauGABA_supVIP_to_semilunar   =20.d0 
      real*8, parameter :: tauGABA_supVIP_to_multipolar=20.d0 
      real*8, parameter :: tauGABA_supVIP_to_deepbask    =20.d0 
      real*8, parameter :: tauGABA_supVIP_to_deepLTS    =20.d0 
      real*8, parameter :: tauGABA_supVIP_to_supVIP      =20.d0 
      real*8, parameter :: tauGABA_supVIP_to_L3pyr   =20.d0 

      real*8, parameter :: tauGABA_deepng_to_semilunar    =6.d0
      real*8, parameter :: tauGABAB_deepng_to_semilunar   =100.d0
      real*8, parameter :: tauGABA_deepng_to_multipolar   =6.d0
      real*8, parameter :: tauGABAB_deepng_to_multipolar   =100.d0
c BUT integration of multipolar may not include GABA-B
      real*8, parameter :: tauGABA_deepng_to_L2pyr    =6.d0
      real*8, parameter :: tauGABA_deepng_to_L3pyr    =6.d0
      real*8, parameter :: tauGABAB_deepng_to_L2pyr   =100.d0
      real*8, parameter :: tauGABAB_deepng_to_L3pyr   =100.d0
      real*8, parameter :: tauGABA_deepng_to_LOT    =6.d0
      real*8, parameter :: tauGABAB_deepng_to_LOT   =100.d0
      real*8, parameter :: tauGABA_deepng_to_deepng       =3.d0
      real*8, parameter :: tauGABA_deepng_to_deepbask     =3.d0

      real*8, parameter :: tauAMPA_placeholder5_to_L2pyr        =2.d0  
      real*8, parameter :: tauNMDA_placeholder5_to_L2pyr      =130.d0
c     real*8, parameter :: tauAMPA_placeholder5_to_supbask         =1.d0  
      real*8, parameter :: tauAMPA_placeholder5_to_supbask      =0.75d0  
      real*8, parameter :: tauNMDA_placeholder5_to_supbask      =100.d0
      real*8, parameter :: tauAMPA_placeholder5_to_supng        =0.75d0  
      real*8, parameter :: tauNMDA_placeholder5_to_supng        =100.d0
      real*8, parameter :: tauAMPA_placeholder5_to_placeholder1   =1.d0  
      real*8, parameter :: tauNMDA_placeholder5_to_placeholder1 =100.d0 
      real*8, parameter :: tauAMPA_placeholder5_to_placeholder2   =1.d0  
      real*8, parameter :: tauNMDA_placeholder5_to_placeholder2 =100.d0 
      real*8, parameter :: tauAMPA_placeholder5_to_LOT       =2.0d0 
      real*8, parameter :: tauNMDA_placeholder5_to_LOT       =130.d0
      real*8, parameter :: tauAMPA_placeholder5_to_semilunar      =2.d0  
      real*8, parameter :: tauNMDA_placeholder5_to_semilunar    =130.d0
      real*8, parameter :: tauAMPA_placeholder5_to_multipolar   =2.d0  
      real*8, parameter :: tauNMDA_placeholder5_to_multipolar   =130.d0
!     real*8, parameter :: tauAMPA_placeholder5_to_deepbask       =1.d0  
      real*8, parameter :: tauAMPA_placeholder5_to_deepbask     =0.75d0  
      real*8, parameter :: tauNMDA_placeholder5_to_deepbask     =100.d0
      real*8, parameter :: tauAMPA_placeholder5_to_deepng       =0.75d0  
      real*8, parameter :: tauNMDA_placeholder5_to_deepng       =100.d0
!     real*8, parameter :: tauAMPA_placeholder5_to_deepLTS      =1.d0  
      real*8, parameter :: tauAMPA_placeholder5_to_deepLTS     =0.75d0  
      real*8, parameter :: tauNMDA_placeholder5_to_deepLTS      =100.d0
      real*8, parameter :: tauAMPA_placeholder5_to_placeholder6  =2.0d0      
      real*8, parameter :: tauNMDA_placeholder5_to_placeholder6 =150.d0
      real*8, parameter :: tauAMPA_placeholder5_to_L3pyr       =2.0d0     
      real*8, parameter :: tauNMDA_placeholder5_to_L3pyr       =130.d0

      real*8, parameter :: tauGABA1_placeholder6_to_placeholder5=3.30d0 
      real*8, parameter :: tauGABA2_placeholder6_to_placeholder5 =10.d0 
      real*8, parameter :: tauGABA1_placeholder6_to_placeholder6 = 9.d0 
      real*8, parameter :: tauGABA2_placeholder6_to_placeholder6=44.5d0 

      real*8, parameter :: tauAMPA_L3pyr_to_L2pyr  =2.d0  
      real*8, parameter :: tauNMDA_L3pyr_to_L2pyr  =130.d0
      real*8, parameter :: tauAMPA_L3pyr_to_supbask   =.8d0  
      real*8, parameter :: tauNMDA_L3pyr_to_supbask   =100.d0
      real*8, parameter :: tauAMPA_L3pyr_to_placeholder1  =.8d0  
      real*8, parameter :: tauNMDA_L3pyr_to_placeholder1  =100.d0 
      real*8, parameter :: tauAMPA_L3pyr_to_placeholder2  =.8d0  
      real*8, parameter :: tauNMDA_L3pyr_to_placeholder2  =100.d0 
      real*8, parameter :: tauAMPA_L3pyr_to_placeholder3  =.8d0  
      real*8, parameter :: tauNMDA_L3pyr_to_placeholder3  =100.d0 
      real*8, parameter :: tauAMPA_L3pyr_to_supLTS    =1.0d0 
      real*8, parameter :: tauNMDA_L3pyr_to_supLTS    =100.d0
      real*8, parameter :: tauAMPA_L3pyr_to_LOT =2.d0  
      real*8, parameter :: tauNMDA_L3pyr_to_LOT =130.d0
      real*8, parameter :: tauAMPA_L3pyr_to_semilunar    =2.d0  
      real*8, parameter :: tauNMDA_L3pyr_to_semilunar    =130.d0
      real*8, parameter :: tauAMPA_L3pyr_to_multipolar    =2.d0  
      real*8, parameter :: tauNMDA_L3pyr_to_multipolar   =130.d0
      real*8, parameter :: tauAMPA_L3pyr_to_deepbask  =.8d0  
      real*8, parameter :: tauNMDA_L3pyr_to_deepbask  =100.d0
      real*8, parameter :: tauAMPA_L3pyr_to_deepng    =.8d0  
      real*8, parameter :: tauNMDA_L3pyr_to_deepng    =100.d0
      real*8, parameter :: tauAMPA_L3pyr_to_deepLTS  =.8d0   
      real*8, parameter :: tauNMDA_L3pyr_to_deepLTS  =100.d0
      real*8, parameter :: tauAMPA_L3pyr_to_supVIP    =1.d0  
      real*8, parameter :: tauNMDA_L3pyr_to_supVIP    =100.d0
      real*8, parameter :: tauAMPA_L3pyr_to_placeholder5  =2.d0  
      real*8, parameter :: tauNMDA_L3pyr_to_placeholder5  =130.d0 
      real*8, parameter :: tauAMPA_L3pyr_to_placeholder6  =2.0d0 
      real*8, parameter :: tauNMDA_L3pyr_to_placeholder6  =100.d0 
      real*8, parameter :: tauAMPA_L3pyr_to_L3pyr =2.d0  
      real*8, parameter :: tauNMDA_L3pyr_to_L3pyr =130.d0 
c End definition of synaptic time constants

c Synaptic conductance scaling factors.
c     real*8 gAMPA_L2pyr_to_L2pyr /15.00d-3/
      real*8 gAMPA_L2pyr_to_L2pyr / 4.000d-3/
c     real*8 gAMPA_L2pyr_to_L2pyr /0.30d-3/
c     real*8 gAMPA_L2pyr_to_L2pyr /0.00d-3/
      real*8 gNMDA_L2pyr_to_L2pyr /0.050d-3/
c     real*8 gNMDA_L2pyr_to_L2pyr /0.000d-3/
c     real*8 gNMDA_L2pyr_to_L2pyr /0.000d-3/
      real*8 gAMPA_L2pyr_to_supbask  /0.00d-3/
      real*8 gNMDA_L2pyr_to_supbask  /0.00d-3/
      real*8 gAMPA_L2pyr_to_supng    /0.30d-3/
      real*8 gNMDA_L2pyr_to_supng    /0.01d-3/
      real*8 gAMPA_L2pyr_to_placeholder1  /0.0d-3/
      real*8 gNMDA_L2pyr_to_placeholder1  /0.00d-3/
      real*8 gAMPA_L2pyr_to_placeholder2  /0.0d-3/
      real*8 gNMDA_L2pyr_to_placeholder2  /0.00d-3/
      real*8 gAMPA_L2pyr_to_placeholder3   /0.0d-3/
      real*8 gNMDA_L2pyr_to_placeholder3   /0.00d-3/
      real*8 gAMPA_L2pyr_to_LOT/0.00d-3/
      real*8 gNMDA_L2pyr_to_LOT/0.00d-3/
      real*8 gAMPA_L2pyr_to_semilunar   /0.20d-3/
      real*8 gNMDA_L2pyr_to_semilunar   /0.01d-3/
      real*8 gAMPA_L2pyr_to_multipolar     /2.00d-3/
      real*8 gNMDA_L2pyr_to_multipolar     /0.01d-3/
      real*8 gAMPA_L2pyr_to_deepbask /2.00d-3/
      real*8 gNMDA_L2pyr_to_deepbask /0.05d-3/
      real*8 gAMPA_L2pyr_to_deepLTS /2.00d-3/
      real*8 gNMDA_L2pyr_to_deepLTS /0.05d-3/
      real*8 gAMPA_L2pyr_to_deepng  /2.00d-3/
      real*8 gNMDA_L2pyr_to_deepng  /0.05d-3/
      real*8 gAMPA_L2pyr_to_supVIP   /0.50d-3/
      real*8 gNMDA_L2pyr_to_supVIP   /0.01d-3/
c     real*8 gAMPA_L2pyr_to_L3pyr /00.00d-3/
      real*8 gAMPA_L2pyr_to_L3pyr /6.00d-3/
      real*8 gNMDA_L2pyr_to_L3pyr /0.05d-3/

      real*8 gGABA_placeholder1_to_L2pyr   /0.0d-3/
      real*8 gGABA_placeholder1_to_placeholder1 /0.0d-3/
      real*8 gGABA_placeholder1_to_supng      /0.0d-3/
      real*8 gGABA_placeholder1_to_placeholder2  /0.0d-3/
      real*8 gGABA_placeholder1_to_placeholder3  /0.0d-3/
      real*8 gGABA_placeholder1_to_LOT  /0.0d-3/

      real*8 gGABA_supng_to_L2pyr     /1.0d-3/
c     real*8 gGABA_supng_to_L2pyr     /0.0d-3/
      real*8 gGABA_supng_to_L3pyr    /1.0d-3/
c     real*8 gGABA_supng_to_L3pyr    /0.0d-3/
      real*8 gGABA_supng_to_semilunar    /1.0d-3/
c     real*8 gGABA_supng_to_semilunar    /0.0d-3/
      real*8 gGABA_supng_to_multipolar       /0.0d-3/
      real*8 gGABA_supng_to_supng        /0.1d-3/
      real*8 gGABA_supng_to_placeholder1      /0.0d-3/

! THESE GABA-B WILL HAVE REL. FAST KINETICS, VS SLOW FROM SUPVIP INTERNEURONS
      real*8 gGABAB_supng_to_L2pyr    /0.010d-3/
      real*8 gGABAB_supng_to_L3pyr   /0.010d-3/
      real*8 gGABAB_supng_to_semilunar     /0.010d-3/
      real*8 gGABAB_supng_to_multipolar  /0.000d-3/

      real*8 gGABAB_supVIP_to_L2pyr    /0.010d-3/
      real*8 gGABAB_supVIP_to_L3pyr   /0.010d-3/
      real*8 gGABAB_supVIP_to_semilunar   /0.010d-3/
      real*8 gGABAB_supVIP_to_multipolar    /0.00d-3/

      real*8 gGABA_placeholder2_to_L2pyr   /0.0d-3/
      real*8 gGABA_placeholder2_to_LOT  /0.0d-3/ 
      real*8 gGABA_placeholder2_to_semilunar     /0.0d-3/
      real*8 gGABA_placeholder2_to_multipolar  /0.0d-3/
      real*8 gGABA_placeholder2_to_L3pyr  /0.0d-3/

      real*8 gGABA_placeholder3_to_L2pyr    /.00d-3/
      real*8 gGABA_placeholder3_to_placeholder1     /.00d-3/
      real*8 gGABA_placeholder3_to_placeholder2     /.00d-3/
      real*8 gGABA_placeholder3_to_placeholder3      /.00d-3/
      real*8 gGABA_placeholder3_to_LOT   /.00d-3/
      real*8 gGABA_placeholder3_to_semilunar      /.00d-3/
      real*8 gGABA_placeholder3_to_multipolar      /.00d-3/
      real*8 gGABA_placeholder3_to_deepbask    /.00d-3/
      real*8 gGABA_placeholder3_to_deepLTS    /.00d-3/
      real*8 gGABA_placeholder3_to_supVIP      /.00d-3/
      real*8 gGABA_placeholder3_to_L3pyr   /.00d-3/

      real*8 gAMPA_LOT_to_L2pyr / 9.0d-3/
c     real*8 gAMPA_LOT_to_L2pyr /00.0d-3/
c     real*8 gNMDA_LOT_to_L2pyr /0.00d-3/
      real*8 gNMDA_LOT_to_L2pyr /0.01d-3/ ! maybe should be large to get
c     real*8 gNMDA_LOT_to_L2pyr /0.00d-3/ ! maybe should be large to get
c        NMDA spikes?
      real*8 gAMPA_LOT_to_placeholder1  /0.0d-3/
      real*8 gNMDA_LOT_to_placeholder1  /.00d-3/
      real*8 gAMPA_LOT_to_placeholder2  /0.0d-3/
      real*8 gNMDA_LOT_to_placeholder2  /.00d-3/
      real*8 gAMPA_LOT_to_placeholder3   /0.0d-3/
      real*8 gNMDA_LOT_to_placeholder3   /.00d-3/
      real*8 gAMPA_LOT_to_LOT/0.0d-3/
      real*8 gNMDA_LOT_to_LOT/0.00d-3/
      real*8 gAMPA_LOT_to_semilunar   /15.0d-3/ ! should be larger than
c  to L2pyr; see Suzuki and Bekkers, J Neurosci 2011
      real*8 gNMDA_LOT_to_semilunar   /0.01d-3/ ! maybe small, as
c apparently no NMDA spikes in SL
      real*8 gAMPA_LOT_to_multipolar   /0.0d-3/
      real*8 gNMDA_LOT_to_multipolar   /0.00d-3/
      real*8 gAMPA_LOT_to_deepbask /0.0d-3/
      real*8 gNMDA_LOT_to_deepbask /.00d-3/
      real*8 gAMPA_LOT_to_deepng   /0.0d-3/
      real*8 gNMDA_LOT_to_deepng   /.00d-3/
      real*8 gAMPA_LOT_to_deepLTS /0.0d-3/
      real*8 gNMDA_LOT_to_deepLTS /.00d-3/
      real*8 gAMPA_LOT_to_supVIP   /4.5d-3/
      real*8 gNMDA_LOT_to_supVIP   /.00d-3/
      real*8 gAMPA_LOT_to_supng    /0.0d-3/
      real*8 gNMDA_LOT_to_supng    /.00d-3/
      real*8 gAMPA_LOT_to_L3pyr/12.0d-3/
      real*8 gNMDA_LOT_to_L3pyr/0.05d-3/
c     real*8 gNMDA_LOT_to_L3pyr/0.00d-3/

      real*8 gAMPA_semilunar_to_L2pyr    /2.0d-3/
c     real*8 gAMPA_semilunar_to_L2pyr    /0.0d-3/
      real*8 gNMDA_semilunar_to_L2pyr    /0.05d-3/
c     real*8 gNMDA_semilunar_to_L2pyr    /0.00d-3/
      real*8 gAMPA_semilunar_to_placeholder1     /0.0d-3/
      real*8 gNMDA_semilunar_to_placeholder1     /0.00d-3/
      real*8 gAMPA_semilunar_to_placeholder2   /0.0d-3/
      real*8 gNMDA_semilunar_to_placeholder2   /0.00d-3/
      real*8 gAMPA_semilunar_to_placeholder3      /0.0d-3/
      real*8 gNMDA_semilunar_to_placeholder3      /0.00d-3/
      real*8 gAMPA_semilunar_to_LOT   /0.0d-3/
      real*8 gNMDA_semilunar_to_LOT   /0.00d-3/
      real*8 gAMPA_semilunar_to_semilunar    /0.1d-3/
      real*8 gNMDA_semilunar_to_semilunar    /0.01d-3/
      real*8 gAMPA_semilunar_to_multipolar    /0.0d-3/
      real*8 gNMDA_semilunar_to_multipolar    /0.00d-3/
      real*8 gAMPA_semilunar_to_deepbask    /1.0d-3/
      real*8 gNMDA_semilunar_to_deepbask    /0.01d-3/
      real*8 gAMPA_semilunar_to_deepng      /1.0d-3/
      real*8 gNMDA_semilunar_to_deepng      /0.01d-3/
      real*8 gAMPA_semilunar_to_deepLTS    /1.0d-3/
      real*8 gNMDA_semilunar_to_deepLTS    /0.01d-3/
      real*8 gAMPA_semilunar_to_supVIP      /0.5d-3/
      real*8 gNMDA_semilunar_to_supVIP      /0.01d-3/
      real*8 gAMPA_semilunar_to_L3pyr   /2.00d-3/
      real*8 gNMDA_semilunar_to_L3pyr   /2.00d-3/

      real*8 gAMPA_multipolar_to_L2pyr    /1.0d-3/
      real*8 gNMDA_multipolar_to_L2pyr    /0.02d-3/
      real*8 gAMPA_multipolar_to_placeholder1     /0.0d-3/
      real*8 gNMDA_multipolar_to_placeholder1     /0.00d-3/
      real*8 gAMPA_multipolar_to_placeholder2   /0.0d-3/
      real*8 gNMDA_multipolar_to_placeholder2   /0.00d-3/
      real*8 gAMPA_multipolar_to_placeholder3      /0.0d-3/
      real*8 gNMDA_multipolar_to_placeholder3      /0.00d-3/
      real*8 gAMPA_multipolar_to_LOT   /0.0d-3/
      real*8 gNMDA_multipolar_to_LOT   /0.00d-3/
      real*8 gAMPA_multipolar_to_semilunar      /0.0d-3/
      real*8 gNMDA_multipolar_to_semilunar      /0.00d-3/
      real*8 gAMPA_multipolar_to_multipolar     /0.60d-3/
      real*8 gNMDA_multipolar_to_multipolar     /0.01d-3/
      real*8 gAMPA_multipolar_to_deepbask    /0.5d-3/
      real*8 gNMDA_multipolar_to_deepbask    /0.05d-3/
      real*8 gAMPA_multipolar_to_deepng      /0.5d-3/
      real*8 gNMDA_multipolar_to_deepng      /0.01d-3/
      real*8 gAMPA_multipolar_to_deepLTS    /0.0d-3/
      real*8 gNMDA_multipolar_to_deepLTS    /0.00d-3/
      real*8 gAMPA_multipolar_to_supVIP      /0.0d-3/
      real*8 gNMDA_multipolar_to_supVIP      /0.00d-3/
      real*8 gAMPA_multipolar_to_L3pyr   /1.0d-3/
      real*8 gNMDA_multipolar_to_L3pyr   /0.05d-3/

      real*8 gGABA_deepbask_to_LOT /0.00d-3/ 
      real*8 gGABA_deepbask_to_semilunar    /1.0d-3/
c     real*8 gGABA_deepbask_to_semilunar    /0.0d-3/
      real*8 gGABA_deepbask_to_multipolar   /1.0d-3/
      real*8 gGABA_deepbask_to_deepbask  /0.5d-3/
      real*8 gGABA_deepbask_to_deepng    /0.2d-3/
      real*8 gGABA_deepbask_to_deepLTS  /0.2d-3/
      real*8 gGABA_deepbask_to_supVIP    /0.0d-3/
      real*8 gGABA_deepbask_to_L2pyr /2.0d-3/
c     real*8 gGABA_deepbask_to_L2pyr /0.0d-3/
      real*8 gGABA_deepbask_to_L3pyr /2.0d-3/
c     real*8 gGABA_deepbask_to_L3pyr /0.0d-3/

      real*8 gGABA_deepng_to_semilunar      /0.0d-3/
c     real*8 gGABA_deepng_to_semilunar      /0.0d-3/
      real*8 gGABA_deepng_to_multipolar        /0.5d-3/
      real*8 gGABA_deepng_to_L2pyr   /0.3d-3/
c     real*8 gGABA_deepng_to_L2pyr   /0.0d-3/
      real*8 gGABA_deepng_to_L3pyr   /0.3d-3/
c     real*8 gGABA_deepng_to_L3pyr   /0.0d-3/
c     real*8 gGABA_deepng_to_L3pyr   /0.0d-3/
      real*8 gGABA_deepng_to_LOT   /0.0d-3/
      real*8 gGABA_deepng_to_deepng      /0.1d-3/
      real*8 gGABA_deepng_to_deepbask    /0.1d-3/

      real*8 gGABAB_deepng_to_semilunar      /0.002d-3/
      real*8 gGABAB_deepng_to_multipolar     /0.002d-3/
      real*8 gGABAB_deepng_to_L2pyr   /0.0020d-3/
      real*8 gGABAB_deepng_to_L3pyr   /0.0020d-3/
      real*8 gGABAB_deepng_to_LOT   /0.000d-3/

      real*8 gGABA_deepLTS_to_L2pyr   /0.1d-3/
c     real*8 gGABA_deepLTS_to_L2pyr   /0.00d-3/
      real*8 gGABA_deepLTS_to_LOT  /0.0d-3/ 
      real*8 gGABA_deepLTS_to_semilunar     /0.1d-3/
c     real*8 gGABA_deepLTS_to_semilunar     /0.0d-3/
      real*8 gGABA_deepLTS_to_multipolar     /0.0d-3/
      real*8 gGABA_deepLTS_to_L3pyr  /0.10d-3/
c     real*8 gGABA_deepLTS_to_L3pyr  /0.0d-3/

c     real*8 gGABA_supVIP_to_L2pyr    /0.00d-3/
      real*8 gGABA_supVIP_to_L2pyr    /1.00d-3/
c     real*8 gGABA_supVIP_to_L2pyr    /.00d-3/
      real*8 gGABA_supVIP_to_placeholder1     /.00d-3/
      real*8 gGABA_supVIP_to_placeholder2     /.00d-3/
      real*8 gGABA_supVIP_to_placeholder3      /.00d-3/
      real*8 gGABA_supVIP_to_LOT   /.00d-3/
c     real*8 gGABA_supVIP_to_semilunar      /0.00d-3/
      real*8 gGABA_supVIP_to_semilunar      /1.00d-3/
c     real*8 gGABA_supVIP_to_semilunar      /.00d-3/
      real*8 gGABA_supVIP_to_multipolar      /.00d-3/
      real*8 gGABA_supVIP_to_deepbask    /.00d-3/
      real*8 gGABA_supVIP_to_deepLTS    /.00d-3/
      real*8 gGABA_supVIP_to_supVIP      /.01d-3/
      real*8 gGABA_supVIP_to_supng       /.01d-3/
c     real*8 gGABA_supVIP_to_L3pyr   /0.00d-3/
      real*8 gGABA_supVIP_to_L3pyr   /1.00d-3/
c     real*8 gGABA_supVIP_to_L3pyr   /.00d-3/

      real*8 gAMPA_placeholder5_to_L2pyr        /0.0d-3/
      real*8 gNMDA_placeholder5_to_L2pyr        /0.00d-3/
      real*8 gAMPA_placeholder5_to_placeholder1      /0.0d-3/
      real*8 gNMDA_placeholder5_to_placeholder1      /.00d-3/
      real*8 gAMPA_placeholder5_to_supng        /0.0d-3/
      real*8 gNMDA_placeholder5_to_supng        /0.0d-3/
      real*8 gAMPA_placeholder5_to_placeholder2 /0.0d-3/
      real*8 gNMDA_placeholder5_to_placeholder2   /.00d-3/
      real*8 gAMPA_placeholder5_to_LOT       /0.0d-3/
      real*8 gNMDA_placeholder5_to_LOT       /.00d-3/
      real*8 gAMPA_placeholder5_to_semilunar      /0.0d-3/
      real*8 gNMDA_placeholder5_to_semilunar      /.00d-3/
      real*8 gAMPA_placeholder5_to_multipolar     /0.0d-3/
      real*8 gNMDA_placeholder5_to_multipolar     /.00d-3/
      real*8 gAMPA_placeholder5_to_deepbask        /0.0d-3/
      real*8 gNMDA_placeholder5_to_deepbask        /.00d-3/
      real*8 gAMPA_placeholder5_to_deepng          /0.0d-3/
      real*8 gNMDA_placeholder5_to_deepng          /0.0d-3/
      real*8 gAMPA_placeholder5_to_deepLTS        /0.0d-3/
      real*8 gNMDA_placeholder5_to_deepLTS        /.00d-3/
      real*8 gAMPA_placeholder5_to_placeholder6     /0.00d-3/   
      real*8 gNMDA_placeholder5_to_placeholder6    /.00d-3/
      real*8 gAMPA_placeholder5_to_L3pyr       /0.0d-3/    
      real*8 gNMDA_placeholder5_to_L3pyr       /.00d-3/

      real*8 gGABAB_placeholder6_to_placeholder5     /0.000d-3/
      real*8 gGABA_placeholder6_to_placeholder5 /0.000d-3/
      real*8 gGABA_placeholder6_to_placeholder6        /0.00d-3/
      real*8 gGABAB_placeholder6_to_placeholder6       /0.000d-3/

      real*8 gAMPA_L3pyr_to_L2pyr  /4.0d-3/
c     real*8 gAMPA_L3pyr_to_L2pyr  /0.0d-3/
c     real*8 gAMPA_L3pyr_to_L2pyr  /15.0d-3/
      real*8 gNMDA_L3pyr_to_L2pyr  /0.05d-3/
c     real*8 gNMDA_L3pyr_to_L2pyr  /0.00d-3/
      real*8 gAMPA_L3pyr_to_placeholder1   /0.0d-3/
      real*8 gNMDA_L3pyr_to_placeholder1   /0.00d-3/
      real*8 gAMPA_L3pyr_to_placeholder2   /0.0d-3/
      real*8 gNMDA_L3pyr_to_placeholder2   /0.00d-3/
      real*8 gAMPA_L3pyr_to_placeholder3    /0.0d-3/
      real*8 gNMDA_L3pyr_to_placeholder3    /0.00d-3/
      real*8 gAMPA_L3pyr_to_LOT /0.0d-3/
      real*8 gNMDA_L3pyr_to_LOT /0.00d-3/
      real*8 gAMPA_L3pyr_to_semilunar    /0.1d-3/
      real*8 gNMDA_L3pyr_to_semilunar    /0.00d-3/
      real*8 gAMPA_L3pyr_to_multipolar      /2.5d-3/
      real*8 gNMDA_L3pyr_to_multipolar      /0.01d-3/
      real*8 gAMPA_L3pyr_to_deepbask  /2.0d-3/
      real*8 gNMDA_L3pyr_to_deepbask  /.01d-3/
      real*8 gAMPA_L3pyr_to_deepng    /1.0d-3/
      real*8 gNMDA_L3pyr_to_deepng    /.01d-3/
      real*8 gAMPA_L3pyr_to_deepLTS  /2.0d-3/
      real*8 gNMDA_L3pyr_to_deepLTS  /.01d-3/
      real*8 gAMPA_L3pyr_to_supVIP    /0.0d-3/
      real*8 gNMDA_L3pyr_to_supVIP    /.00d-3/
      real*8 gAMPA_L3pyr_to_placeholder5       /.00d-3/ 
      real*8 gNMDA_L3pyr_to_placeholder5       /.000d-3/
      real*8 gAMPA_L3pyr_to_placeholder6       /0.0d-3/
      real*8 gNMDA_L3pyr_to_placeholder6       /0.00d-3/
      real*8 gAMPA_L3pyr_to_L3pyr /4.0d-3/
c     real*8 gAMPA_L3pyr_to_L3pyr /00.0d-3/
      real*8 gNMDA_L3pyr_to_L3pyr /0.01d-3/
c End defining synaptic conductance scaling factors

c Begin definition of compartments where synaptic connections
c can form.
       INTEGER compallow_L2pyr_to_L2pyr 
     &  (ncompallow_L2pyr_to_L2pyr)
     &  /2,3,4,5,6,7,8,9,10,11,12,13,22,23,24,25,
     &  14,15,16,17,18,19,20,21,
     &  26,33,34,35,36,37,39,40,41,42,43,44/
       INTEGER compallow_L2pyr_to_placeholder1  
     &  (ncompallow_L2pyr_to_placeholder1  )
     &  /1/
       INTEGER compallow_L2pyr_to_supng    
     &  (ncompallow_L2pyr_to_supng    )
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_L2pyr_to_placeholder2  
     &  (ncompallow_L2pyr_to_placeholder2  )
     &  /1/
       INTEGER compallow_L2pyr_to_placeholder3   
     &  (ncompallow_L2pyr_to_placeholder3   )
     &  /1/
       INTEGER compallow_L2pyr_to_LOT
     &  (ncompallow_L2pyr_to_LOT)
     &  /1/
       INTEGER compallow_L2pyr_to_semilunar   
     &  (ncompallow_L2pyr_to_semilunar   )
     &  /38,39,40,41,42,43,44/
       INTEGER compallow_L2pyr_to_multipolar     
     &  (ncompallow_L2pyr_to_multipolar     )
     &  /3,4,16,17,29,30,42,43,5,7,18,19,31,32,44,45,
     &  6,20,33,46,21,22,23,24/
       INTEGER compallow_L2pyr_to_deepbask 
     &  (ncompallow_L2pyr_to_deepbask )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_L2pyr_to_deepLTS 
     &  (ncompallow_L2pyr_to_deepLTS )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
        INTEGER compallow_L2pyr_to_deepng   
     &    (ncompallow_L2pyr_to_deepng  ) 
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_L2pyr_to_supVIP   
     &  (ncompallow_L2pyr_to_supVIP   )
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_L2pyr_to_L3pyr
     &  (ncompallow_L2pyr_to_L3pyr)
     &  /2,3,4,5,6,7,8,9,10,11,12,13,22,23,24,25,
     &  14,15,16,17,18,19,20,21,
     &  26,33,34,35,36,37,39,40,41,42,43,44/

       INTEGER compallow_placeholder1_to_L2pyr
     &  (ncompallow_placeholder1_to_L2pyr)
     & /1/
       INTEGER compallow_placeholder1_to_placeholder1  
     &  (ncompallow_placeholder1_to_placeholder1  )
     &  /1/
       INTEGER compallow_placeholder1_to_supng    
     &  (ncompallow_placeholder1_to_supng    )
     & /1/             
       INTEGER compallow_placeholder1_to_placeholder2  
     &  (ncompallow_placeholder1_to_placeholder2  )
     &  /1/
       INTEGER compallow_placeholder1_to_placeholder3   
     &  (ncompallow_placeholder1_to_placeholder3   )
     &  /1/
       INTEGER compallow_placeholder1_to_LOT
     &  (ncompallow_placeholder1_to_LOT)
     &  /1/

       INTEGER compallow_supng_to_L2pyr 
     &  (ncompallow_supng_to_L2pyr )
     & /45,46,47,48,49,50,51,52,53,54,55,56,57,58,
     &  59,60,61,62,63,64,65,66,67,68/
       INTEGER compallow_supng_to_L3pyr
     &  (ncompallow_supng_to_L3pyr)
     & /45,46,47,48,49,50,51,52,53,54,55,56,57,58,
     &  59,60,61,62,63,64,65,66,67,68/
       INTEGER compallow_supng_to_semilunar   
     &  (ncompallow_supng_to_semilunar   )
     & /45,46,47,48,49,50,51,52,53,54,55,56,57,58,
     &  59,60,61,62,63,64,65,66,67,68/
       INTEGER compallow_supng_to_multipolar     
     &  (ncompallow_supng_to_multipolar     )
     & /1/
       INTEGER compallow_supng_to_supng    
     &  (ncompallow_supng_to_supng    )
     & /2,1,28,41/
       INTEGER compallow_supng_to_placeholder1  
     &  (ncompallow_supng_to_placeholder1  )
     & /1/

       INTEGER compallow_placeholder3_to_L2pyr
     &  (ncompallow_placeholder3_to_L2pyr)
     & /1/
       INTEGER compallow_placeholder3_to_placeholder1  
     &  (ncompallow_placeholder3_to_placeholder1)  
     & /1/
       INTEGER compallow_placeholder3_to_placeholder2  
     &  (ncompallow_placeholder3_to_placeholder2)  
     & /1/ 
       INTEGER compallow_placeholder3_to_placeholder3   
     &  (ncompallow_placeholder3_to_placeholder3 )  
     & /1/
       INTEGER compallow_placeholder3_to_LOT
     &  (ncompallow_placeholder3_to_LOT)
     & /1/
       INTEGER compallow_placeholder3_to_semilunar   
     &  (ncompallow_placeholder3_to_semilunar   )
     & / 1/ 
       INTEGER compallow_placeholder3_to_multipolar     
     &  (ncompallow_placeholder3_to_multipolar     )
     & / 1/
       INTEGER compallow_placeholder3_to_deepbask 
     &  (ncompallow_placeholder3_to_deepbask )
     & / 1/
       INTEGER compallow_placeholder3_to_deepLTS 
     &  (ncompallow_placeholder3_to_deepLTS )
     & / 1/
       INTEGER compallow_placeholder3_to_supVIP   
     &  (ncompallow_placeholder3_to_supVIP   )
     & / 1/
       INTEGER compallow_placeholder3_to_L3pyr
     &  (ncompallow_placeholder3_to_L3pyr)
     & / 1/

       INTEGER compallow_LOT_to_L2pyr
     &   (ncompallow_LOT_to_L2pyr)
     & / 45,46,47,48,49,50,51,52,53,54,55,56,57,58,
     &   59,60,61,62,63,64,65,66,67,68/                        
       INTEGER compallow_LOT_to_placeholder1  
     &   (ncompallow_LOT_to_placeholder1  )
     &  /1 /
       INTEGER compallow_LOT_to_placeholder2  
     &   (ncompallow_LOT_to_placeholder2  )
     &  /1/
       INTEGER compallow_LOT_to_placeholder3   
     &   (ncompallow_LOT_to_placeholder3   )
     &  /1/
       INTEGER compallow_LOT_to_LOT
     &   (ncompallow_LOT_to_LOT)
     &  /1/
       INTEGER compallow_LOT_to_semilunar   
     &   (ncompallow_LOT_to_semilunar   )
     & / 45,46,47,48,49,50,51,52,53,54,55,56,57,58,
     &   59,60,61,62,63,64,65,66,67,68/                        
       INTEGER compallow_LOT_to_multipolar     
     &   (ncompallow_LOT_to_multipolar     )
     &  /1/
       INTEGER compallow_LOT_to_deepbask 
     &   (ncompallow_LOT_to_deepbask )
     &  /1/
       INTEGER compallow_LOT_to_deepng   
     &   (ncompallow_LOT_to_deepng   )
     &  /1/
       INTEGER compallow_LOT_to_deepLTS 
     &   (ncompallow_LOT_to_deepLTS )
     &  /1/
       INTEGER compallow_LOT_to_supVIP   
     &   (ncompallow_LOT_to_supVIP   )
c    & / 45,46,47,48,49,50,51,52,53,54,55,56,57,58,
c    &   59,60,61,62,63,64,65,66,67,68/                        
     & /45,46,47,48,49,50,51,52,53/ ! corrected 22Mar2021
       INTEGER compallow_LOT_to_supng    
     &   (ncompallow_LOT_to_supng    )
c    & / 45,46,47,48,49,50,51,52,53,54,55,56,57,58,
c    &   59,60,61,62,63,64,65,66,67,68/                        
     & /45,46,47,48,49,50,51,52,53/
       INTEGER compallow_LOT_to_L3pyr
     &   (ncompallow_LOT_to_L3pyr)
     & / 45,46,47,48,49,50,51,52,53,54,55,56,57,58,
     &   59,60,61,62,63,64,65,66,67,68/                        

       INTEGER compallow_semilunar_to_L2pyr
     &   (ncompallow_semilunar_to_L2pyr)
     & /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
     & 20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,
     & 36,37,38,39,40,41,42,43,44/ 
       INTEGER compallow_semilunar_to_placeholder1  
     &   (ncompallow_semilunar_to_placeholder1)  
     &  /1/
       INTEGER compallow_semilunar_to_placeholder2  
     &   (ncompallow_semilunar_to_placeholder2)  
     &  /1/
       INTEGER compallow_semilunar_to_placeholder3   
     &   (ncompallow_semilunar_to_placeholder3 )  
     &  /1/
       INTEGER compallow_semilunar_to_LOT
     &   (ncompallow_semilunar_to_LOT) 
     &  /1/
       INTEGER compallow_semilunar_to_semilunar   
     &   (ncompallow_semilunar_to_semilunar)    
     &   /10,11,12,13,22,23,24,25,34,35,36,37,
     & 38,39,40/
       INTEGER compallow_semilunar_to_multipolar     
     &   (ncompallow_semilunar_to_multipolar  )    
     &  /1/
       INTEGER compallow_semilunar_to_deepbask 
     &   (ncompallow_semilunar_to_deepbask)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_semilunar_to_deepng   
     &   (ncompallow_semilunar_to_deepng  )  
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/
       INTEGER compallow_semilunar_to_deepLTS 
     &   (ncompallow_semilunar_to_deepLTS)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_semilunar_to_supVIP   
     &   (ncompallow_semilunar_to_supVIP  )  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_semilunar_to_L3pyr
     &   (ncompallow_semilunar_to_L3pyr) 
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &   37,38,39,40,41,42,43,44/

       INTEGER compallow_multipolar_to_L2pyr
     &   (ncompallow_multipolar_to_L2pyr)
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &   37/
       INTEGER compallow_multipolar_to_placeholder1  
     &   (ncompallow_multipolar_to_placeholder1)  
     &  /1/
       INTEGER compallow_multipolar_to_placeholder2  
     &   (ncompallow_multipolar_to_placeholder2)  
     &  /1/
       INTEGER compallow_multipolar_to_placeholder3   
     &   (ncompallow_multipolar_to_placeholder3 )  
     &  /1/
       INTEGER compallow_multipolar_to_LOT
     &   (ncompallow_multipolar_to_LOT) 
     &  /1/
       INTEGER compallow_multipolar_to_semilunar   
     &   (ncompallow_multipolar_to_semilunar)    
     &  /1/
       INTEGER compallow_multipolar_to_multipolar   
     &   (ncompallow_multipolar_to_multipolar)    
     &  /3,4,16,17,29,30,42,43,5,7,18,19,31,32,44,45,
     &  6,20,33,46,21,22,23,24/
       INTEGER compallow_multipolar_to_deepbask 
     &   (ncompallow_multipolar_to_deepbask)  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_multipolar_to_deepng   
     &   (ncompallow_multipolar_to_deepng  )  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_multipolar_to_deepLTS 
     &   (ncompallow_multipolar_to_deepLTS)  
     &  /1/
       INTEGER compallow_multipolar_to_supVIP   
     &   (ncompallow_multipolar_to_supVIP  )  
     &  /1/
       INTEGER compallow_multipolar_to_L3pyr
     &   (ncompallow_multipolar_to_L3pyr) 
     &  /2,3,4,5,6,7,8,9,10,11,12,13,22,23,24,25,
     &  14,15,16,17,18,19,20,21,
     &  26,33,34,35,36,37,39,40,41,42,43,44/

       INTEGER compallow_deepbask_to_LOT
     &   (ncompallow_deepbask_to_LOT)
     &  /1/
       INTEGER compallow_deepbask_to_semilunar   
     &   (ncompallow_deepbask_to_semilunar)   
     & / 1,2,3,4,5,6,7,8,9,38,39/
       INTEGER compallow_deepbask_to_multipolar   
     &   (ncompallow_deepbask_to_multipolar)   
     & / 1,2,15,28,41,3,4,16,17,29,30,42,43,5,6,
     & 18,19,31,32,44,45,7,8,9/
       INTEGER compallow_deepbask_to_deepbask 
     &   (ncompallow_deepbask_to_deepbask) 
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_deepbask_to_deepng   
     &   (ncompallow_deepbask_to_deepng  ) 
     &  /2,15,28,41/
       INTEGER compallow_deepbask_to_deepLTS 
     &   (ncompallow_deepbask_to_deepLTS) 
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_deepbask_to_supVIP   
     &   (ncompallow_deepbask_to_supVIP )  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
       INTEGER compallow_deepbask_to_L2pyr
     &   (ncompallow_deepbask_to_L2pyr)
     &  /1,2,3,4,5,6,7,8,9,38,39/
       INTEGER compallow_deepbask_to_L3pyr
     &   (ncompallow_deepbask_to_L3pyr)
     &  /1,2,3,4,5,6,7,8,9,38,39/

       INTEGER compallow_deepng_to_semilunar    
     &  (ncompallow_deepng_to_semilunar   )
     & /1,10,11,12,13,22,23,24,25,34,35,36,37,38,39,40,
     & 41,42,43,44,45,46,47,48,49,50,51,52,53,54,56,57,58/
       INTEGER compallow_deepng_to_multipolar    
     &  (ncompallow_deepng_to_multipolar   )
     & / 1,2,15,28,41,3,4,16,17,29,30,42,43,5,6,
     & 18,19,31,32,44,45,7,8,9/
       INTEGER compallow_deepng_to_L2pyr 
     &  (ncompallow_deepng_to_L2pyr)
     & /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34/
       INTEGER compallow_deepng_to_L3pyr 
     &  (ncompallow_deepng_to_L3pyr)
     & /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34/
       INTEGER compallow_deepng_to_LOT 
     &  (ncompallow_deepng_to_LOT)
     &  /1/
       INTEGER compallow_deepng_to_deepng    
     &  (ncompallow_deepng_to_deepng   )
     &  /2,15,28,41/
       INTEGER compallow_deepng_to_deepbask  
     &  (ncompallow_deepng_to_deepbask )
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/

        INTEGER compallow_deepLTS_to_L2pyr
     &   (ncompallow_deepLTS_to_L2pyr)
     & /45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
     &  61,62,63,64,65,66,67,68/
        INTEGER compallow_deepLTS_to_LOT
     &   (ncompallow_deepLTS_to_LOT) /1/
        INTEGER compallow_deepLTS_to_semilunar
     &   (ncompallow_deepLTS_to_semilunar) 
     & /45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
     &  61,62,63,64,65,66,67,68/
        INTEGER compallow_deepLTS_to_multipolar  
     &    (ncompallow_deepLTS_to_multipolar) /1/
        INTEGER compallow_deepLTS_to_L3pyr
     &    (ncompallow_deepLTS_to_L3pyr)
     & /45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
     &  61,62,63,64,65,66,67,68/

       INTEGER compallow_supVIP_to_L2pyr
     &   (ncompallow_supVIP_to_L2pyr)
     & /45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
     &  61,62,63,64,65,66,67,68/
       INTEGER compallow_supVIP_to_placeholder1  
     &   (ncompallow_supVIP_to_placeholder1)  
     & /1/
       INTEGER compallow_supVIP_to_placeholder2  
     &   (ncompallow_supVIP_to_placeholder2)  
     & /1/
       INTEGER compallow_supVIP_to_placeholder3   
     &   (ncompallow_supVIP_to_placeholder3)   
     & /1/
       INTEGER compallow_supVIP_to_supng    
     &   (ncompallow_supVIP_to_supng)   
     & / 8,9,10,11,12,21,22,23,24,25,34,35,36,37,38,
     &   47,48,49,50,51/ 
       INTEGER compallow_supVIP_to_LOT
     &   (ncompallow_supVIP_to_LOT)
     & /1/
       INTEGER compallow_supVIP_to_semilunar   
     &   (ncompallow_supVIP_to_semilunar)    
     & /45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
     &  61,62,63,64,65,66,67,68/
       INTEGER compallow_supVIP_to_multipolar   
     &   (ncompallow_supVIP_to_multipolar)    
     & /1/
       INTEGER compallow_supVIP_to_deepbask 
     &   (ncompallow_supVIP_to_deepbask)  
     & /2,15,28,41/
       INTEGER compallow_supVIP_to_deepLTS 
     &   (ncompallow_supVIP_to_deepLTS)  
     & /2,15,28,41/
       INTEGER compallow_supVIP_to_supVIP  
     &   (ncompallow_supVIP_to_supVIP)   
     & /2,15,28,41/
       INTEGER compallow_supVIP_to_L3pyr
     &   (ncompallow_supVIP_to_L3pyr) 
     & /45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
     &  61,62,63,64,65,66,67,68/

       INTEGER compallow_placeholder5_to_L2pyr
     &   (ncompallow_placeholder5_to_L2pyr)
     &  /1/
       INTEGER compallow_placeholder5_to_placeholder1  
     &   (ncompallow_placeholder5_to_placeholder1)  
     &  /1/
       INTEGER compallow_placeholder5_to_supng    
     &   (ncompallow_placeholder5_to_supng  )  
     &  /1/
       INTEGER compallow_placeholder5_to_placeholder2  
     &   (ncompallow_placeholder5_to_placeholder2)  
     &  /1/
       INTEGER compallow_placeholder5_to_LOT
     &   (ncompallow_placeholder5_to_LOT)
     &  /1/
       INTEGER compallow_placeholder5_to_semilunar   
     &   (ncompallow_placeholder5_to_semilunar)   
     &  /1/
       INTEGER compallow_placeholder5_to_multipolar     
     &   (ncompallow_placeholder5_to_multipolar  )   
     &  /1/
       INTEGER compallow_placeholder5_to_deepbask 
     &   (ncompallow_placeholder5_to_deepbask) 
     &  /1/
       INTEGER compallow_placeholder5_to_deepng   
     &   (ncompallow_placeholder5_to_deepng  ) 
     &  /1/
       INTEGER compallow_placeholder5_to_deepLTS 
     &   (ncompallow_placeholder5_to_deepLTS) 
     &  /1/
       INTEGER compallow_placeholder5_to_placeholder6      
     &   (ncompallow_placeholder5_to_placeholder6)      
     &  /1/
       INTEGER compallow_placeholder5_to_L3pyr
     &   (ncompallow_placeholder5_to_L3pyr)
     &  /1/

       INTEGER compallow_placeholder6_to_placeholder5
     &   (ncompallow_placeholder6_to_placeholder5)
     &  /1/
       INTEGER compallow_placeholder6_to_placeholder6
     &   (ncompallow_placeholder6_to_placeholder6)
     &  /1/

        INTEGER compallow_L3pyr_to_L2pyr
     &    (ncompallow_L3pyr_to_L2pyr)
     &   /2,3,4,5,6,7,8,9,14,15,16,17,18,19,20,21,
     & 39,40,41,42,43,44,10,11,12,13,22,23,24,25,26,
     & 33,34,35,36,37/ 
        INTEGER compallow_L3pyr_to_placeholder1  
     &    (ncompallow_L3pyr_to_placeholder1)  
     &  /1/
        INTEGER compallow_L3pyr_to_placeholder2  
     &    (ncompallow_L3pyr_to_placeholder2)  
     &  /1/
        INTEGER compallow_L3pyr_to_placeholder3   
     &    (ncompallow_L3pyr_to_placeholder3)   
     &  /1/
        INTEGER compallow_L3pyr_to_LOT
     &    (ncompallow_L3pyr_to_LOT)
     &  /1/
        INTEGER compallow_L3pyr_to_semilunar   
     &    (ncompallow_L3pyr_to_semilunar)   
     &  /38,39,40,41,42,43,44/
        INTEGER compallow_L3pyr_to_multipolar   
     &    (ncompallow_L3pyr_to_multipolar)   
     &  /3,4,16,17,29,30,42,43,5,7,18,19,31,32,44,45,
     &  6,20,33,46,21,22,23,24/
        INTEGER compallow_L3pyr_to_deepbask 
     &    (ncompallow_L3pyr_to_deepbask) 
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
        INTEGER compallow_L3pyr_to_deepng   
     &    (ncompallow_L3pyr_to_deepng  ) 
     &  /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     &  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &  37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53/
        INTEGER compallow_L3pyr_to_deepLTS 
     &    (ncompallow_L3pyr_to_deepLTS) 
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
        INTEGER compallow_L3pyr_to_supVIP   
     &    (ncompallow_L3pyr_to_supVIP  )  
     &  /5,6,7,8,9,10,18,19,20,21,22,23,31,32,33,34,35,36,
     &   44,45,46,47,48,49/
        INTEGER compallow_L3pyr_to_placeholder5      
     &    (ncompallow_L3pyr_to_placeholder5)      
     &  /1/
        INTEGER compallow_L3pyr_to_placeholder6      
     &    (ncompallow_L3pyr_to_placeholder6)      
     & /1/
        INTEGER compallow_L3pyr_to_L3pyr
     &    (ncompallow_L3pyr_to_L3pyr)
     &   /2,3,4,5,6,7,8,9,14,15,16,17,18,19,20,21,
     & 39,40,41,42,43,44,10,11,12,13,22,23,24,25,26,
     & 33,34,35,36,37/ 


c Maps of synaptic connectivity.  For simplicity, all contacts
c only made to one compartment.  Axoaxonic cells forced to contact 
c initial axon segments; other compartments will be listed in arrays.
        INTEGER 
     & map_L2pyr_to_L2pyr(num_L2pyr_to_L2pyr,
     &                           num_L2pyr),
     & map_L2pyr_to_placeholder1(num_L2pyr_to_placeholder1,  
     &                           num_placeholder1), 
     & map_L2pyr_to_supng  (num_L2pyr_to_supng  ,  
     &                           num_supng  ), 
     & map_L2pyr_to_placeholder2(num_L2pyr_to_placeholder2, 
     &                           num_placeholder2),
     & map_L2pyr_to_placeholder3(num_L2pyr_to_placeholder3,   
     &                           num_placeholder3),
     & map_L2pyr_to_LOT(num_L2pyr_to_LOT,
     &                           num_LOT),
     & map_L2pyr_to_semilunar(num_L2pyr_to_semilunar,
     &                           num_semilunar),  
     & map_L2pyr_to_multipolar(num_L2pyr_to_multipolar,
     &                           num_multipolar), 
     & map_L2pyr_to_deepbask(num_L2pyr_to_deepbask,
     &                           num_deepbask), 
     & map_L2pyr_to_deepLTS(num_L2pyr_to_deepLTS,
     &                           num_deepLTS), 
     & map_L2pyr_to_deepng  (num_L2pyr_to_deepng  ,
     &                             num_deepng  ), 
     & map_L2pyr_to_supVIP (num_L2pyr_to_supVIP ,
     &                           num_supVIP ), 
     & map_L2pyr_to_L3pyr(num_L2pyr_to_L3pyr,
     &                           num_L3pyr) 
              INTEGER
     & map_placeholder1_to_L2pyr(num_placeholder1_to_L2pyr,
     &                           num_L2pyr),  
     &map_placeholder1_to_placeholder1(num_placeholder1_to_placeholder1,
     &                           num_placeholder1), 
     & map_placeholder1_to_supng  (num_placeholder1_to_supng  ,
     &                           num_supng  ),  
     &map_placeholder1_to_placeholder2(num_placeholder1_to_placeholder2,
     &                           num_placeholder2),
     &map_placeholder1_to_placeholder3(num_placeholder1_to_placeholder3,
     &                           num_placeholder3),  
     & map_placeholder1_to_LOT(num_placeholder1_to_LOT,
     &                           num_LOT)  
              INTEGER
     & map_supng_to_L2pyr  (num_supng_to_L2pyr ,
     &                           num_L2pyr),
     & map_supng_to_L3pyr (num_supng_to_L3pyr,
     &                           num_L3pyr),
     & map_supng_to_semilunar    (num_supng_to_semilunar   ,
     &                           num_semilunar   ),
     & map_supng_to_multipolar    (num_supng_to_multipolar   ,
     &                           num_multipolar   ),
     & map_supng_to_supng     (num_supng_to_supng    ,
     &                           num_supng    ),
     & map_supng_to_placeholder1   (num_supng_to_placeholder1  ,
     &                           num_placeholder1  ), 

     & map_placeholder2_to_L2pyr(num_placeholder2_to_L2pyr,
     &                           num_L2pyr), 
     & map_placeholder2_to_LOT(num_placeholder2_to_LOT,
     &                           num_LOT),
     & map_placeholder2_to_semilunar(num_placeholder2_to_semilunar,
     &                           num_semilunar),  
     &map_placeholder2_to_multipolar(num_placeholder2_to_multipolar,
     &                           num_multipolar), 
     & map_placeholder2_to_L3pyr(num_placeholder2_to_L3pyr,
     &                           num_L3pyr), 
     & map_placeholder3_to_L2pyr(num_placeholder3_to_L2pyr,
     &                           num_L2pyr),  
     &map_placeholder3_to_placeholder1(num_placeholder3_to_placeholder1,
     &                           num_placeholder1),  
     &map_placeholder3_to_placeholder2(num_placeholder3_to_placeholder2,
     &                           num_placeholder2), 
     &map_placeholder3_to_placeholder3(num_placeholder3_to_placeholder3,
     &                           num_placeholder3), 
     & map_placeholder3_to_LOT(num_placeholder3_to_LOT,
     &                           num_LOT), 
     & map_placeholder3_to_semilunar(num_placeholder3_to_semilunar,
     &                           num_semilunar),   
     &map_placeholder3_to_multipolar(num_placeholder3_to_multipolar,
     &                           num_multipolar),  
     & map_placeholder3_to_deepbask(num_placeholder3_to_deepbask,
     &                           num_deepbask), 
     & map_placeholder3_to_deepLTS(num_placeholder3_to_deepLTS,
     &                           num_deepLTS), 
     & map_placeholder3_to_supVIP (num_placeholder3_to_supVIP ,
     &                           num_supVIP ), 
     & map_placeholder3_to_L3pyr(num_placeholder3_to_L3pyr,
     &                           num_L3pyr), 
     & map_LOT_to_L2pyr(num_LOT_to_L2pyr,
     &                           num_L2pyr),
     & map_LOT_to_placeholder1(num_LOT_to_placeholder1,
     &                           num_placeholder1) 
               INTEGER
     & map_LOT_to_placeholder2(num_LOT_to_placeholder2,
     &                           num_placeholder2),
     & map_LOT_to_placeholder3(num_LOT_to_placeholder3,
     &                           num_placeholder3), 
     & map_LOT_to_LOT(num_LOT_to_LOT,
     &                           num_LOT),
     & map_LOT_to_semilunar(num_LOT_to_semilunar,
     &                           num_semilunar),  
     & map_LOT_to_multipolar(num_LOT_to_multipolar,
     &                           num_multipolar), 
     & map_LOT_to_deepbask(num_LOT_to_deepbask,
     &                           num_deepbask), 
     & map_LOT_to_deepng  (num_LOT_to_deepng  ,
     &                           num_deepng  ), 
     & map_LOT_to_deepLTS(num_LOT_to_deepLTS,
     &                           num_deepLTS),
     & map_LOT_to_supVIP (num_LOT_to_supVIP ,
     &                           num_supVIP ),
     & map_LOT_to_supng  (num_LOT_to_supng  ,
     &                           num_supng  ),
     & map_LOT_to_L3pyr(num_LOT_to_L3pyr,
     &                           num_L3pyr),

     & map_semilunar_to_L2pyr(num_semilunar_to_L2pyr,
     &                           num_L2pyr),   
     & map_semilunar_to_placeholder1(num_semilunar_to_placeholder1,
     &                           num_placeholder1),  
     & map_semilunar_to_placeholder2(num_semilunar_to_placeholder2,
     &                           num_placeholder2), 
     & map_semilunar_to_placeholder3(num_semilunar_to_placeholder3,
     &                           num_placeholder3), 
     & map_semilunar_to_LOT(num_semilunar_to_LOT,
     &                           num_LOT), 
     & map_semilunar_to_semilunar(num_semilunar_to_semilunar,
     &                           num_semilunar),   
     & map_semilunar_to_multipolar(num_semilunar_to_multipolar,
     &                           num_multipolar),  
     & map_semilunar_to_deepbask(num_semilunar_to_deepbask,
     &                           num_deepbask), 
     & map_semilunar_to_deepng  (num_semilunar_to_deepng  ,
     &                           num_deepng  ), 
     & map_semilunar_to_deepLTS(num_semilunar_to_deepLTS,
     &                           num_deepLTS),  
     & map_semilunar_to_supVIP (num_semilunar_to_supVIP ,
     &                           num_supVIP ),  
     & map_semilunar_to_L3pyr(num_semilunar_to_L3pyr,
     &                           num_L3pyr), 
     & map_multipolar_to_L2pyr(num_multipolar_to_L2pyr,
     &                           num_L2pyr), 
     &map_multipolar_to_placeholder1(num_multipolar_to_placeholder1,
     &                           num_placeholder1),  
     &map_multipolar_to_placeholder2(num_multipolar_to_placeholder2,
     &                           num_placeholder2),   
     &map_multipolar_to_placeholder3(num_multipolar_to_placeholder3,
     &                           num_placeholder3)     
            INTEGER
     & map_multipolar_to_LOT(num_multipolar_to_LOT,
     &                           num_LOT), 
     & map_multipolar_to_semilunar(num_multipolar_to_semilunar,
     &                           num_semilunar),   
     &map_multipolar_to_multipolar(num_multipolar_to_multipolar,
     &                           num_multipolar),     
     & map_multipolar_to_deepbask(num_multipolar_to_deepbask,
     &                           num_deepbask),  
     & map_multipolar_to_deepng  (num_multipolar_to_deepng  ,
     &                           num_deepng  ),  
     & map_multipolar_to_deepLTS(num_multipolar_to_deepLTS,
     &                           num_deepLTS),   
     & map_multipolar_to_supVIP (num_multipolar_to_supVIP ,
     &                           num_supVIP ),   
     & map_multipolar_to_L3pyr(num_multipolar_to_L3pyr,
     &                           num_L3pyr),  
     & map_deepbask_to_LOT(num_deepbask_to_LOT,
     &                           num_LOT), 
     & map_deepbask_to_semilunar(num_deepbask_to_semilunar,
     &                           num_semilunar),   
     & map_deepbask_to_multipolar(num_deepbask_to_multipolar,
     &                           num_multipolar),  
     & map_deepbask_to_deepbask(num_deepbask_to_deepbask,
     &                           num_deepbask), 
     & map_deepbask_to_deepng  (num_deepbask_to_deepng  ,
     &                           num_deepng  ), 
     & map_deepbask_to_deepLTS(num_deepbask_to_deepLTS,
     &                           num_deepLTS),  
     & map_deepbask_to_supVIP (num_deepbask_to_supVIP ,
     &                           num_supVIP )  
                INTEGER
     & map_deepbask_to_L2pyr(num_deepbask_to_L2pyr,
     &                           num_L3pyr), 
     & map_deepbask_to_L3pyr(num_deepbask_to_L3pyr,
     &                           num_L3pyr), 
     & map_deepng_to_semilunar     (num_deepng_to_semilunar     ,
     &                           num_semilunar      ),
     & map_deepng_to_multipolar     (num_deepng_to_multipolar     ,
     &                           num_multipolar      ),
     & map_deepng_to_L2pyr  (num_deepng_to_L2pyr  ,
     &                           num_L3pyr   ),
     & map_deepng_to_L3pyr  (num_deepng_to_L3pyr  ,
     &                           num_L3pyr   ),
     & map_deepng_to_LOT  (num_deepng_to_LOT  ,
     &                           num_LOT   ),
     & map_deepng_to_deepng     (num_deepng_to_deepng     ,
     &                           num_deepng      ),
     & map_deepng_to_deepbask   (num_deepng_to_deepbask   ,
     &                           num_deepbask    ) 

                INTEGER
     & map_deepLTS_to_L2pyr(num_deepLTS_to_L2pyr,
     &                           num_L2pyr), 
     & map_deepLTS_to_LOT(num_deepLTS_to_LOT,
     &                           num_LOT),
     & map_deepLTS_to_semilunar(num_deepLTS_to_semilunar,
     &                           num_semilunar), 
     & map_deepLTS_to_multipolar(num_deepLTS_to_multipolar,
     &                           num_multipolar),    
     & map_deepLTS_to_L3pyr(num_deepLTS_to_L3pyr,
     &                           num_L3pyr)

                 INTEGER
     & map_supVIP_to_L2pyr(num_supVIP_to_L2pyr,
     &                           num_L2pyr), 
     & map_supVIP_to_placeholder1(num_supVIP_to_placeholder1,
     &                           num_placeholder1),  
     & map_supVIP_to_placeholder2(num_supVIP_to_placeholder2,
     &                           num_placeholder2), 
     & map_supVIP_to_placeholder3(num_supVIP_to_placeholder3,
     &                           num_placeholder3), 
     & map_supVIP_to_supng (num_supVIP_to_supng ,
     &                           num_supng ), 
     & map_supVIP_to_LOT(num_supVIP_to_LOT,
     &                           num_LOT),
     & map_supVIP_to_semilunar(num_supVIP_to_semilunar,
     &                           num_semilunar),  
     & map_supVIP_to_multipolar(num_supVIP_to_multipolar,
     &                            num_multipolar), 
     & map_supVIP_to_deepbask(num_supVIP_to_deepbask,
     &                            num_deepbask), 
     & map_supVIP_to_deepLTS(num_supVIP_to_deepLTS,
     &                            num_deepLTS),  
     & map_supVIP_to_supVIP(num_supVIP_to_supVIP,
     &                            num_supVIP),  
     & map_supVIP_to_L3pyr(num_supVIP_to_L3pyr,
     &                            num_L3pyr), 

     & map_placeholder5_to_L2pyr(num_placeholder5_to_L2pyr,
     &                            num_L2pyr),     
     & map_placeholder5_to_placeholder1(num_placeholder5_to_placeholder1
     &                      ,     num_placeholder1)    
               INTEGER
     & map_placeholder5_to_supng  (num_placeholder5_to_supng  ,
     &                            num_supng  ),   
     &map_placeholder5_to_placeholder2(num_placeholder5_to_placeholder2,
     &  num_placeholder2),
     & map_placeholder5_to_LOT(num_placeholder5_to_LOT,num_LOT),
     & map_placeholder5_to_semilunar(num_placeholder5_to_semilunar,
     &   num_semilunar),
     &map_placeholder5_to_multipolar(num_placeholder5_to_multipolar,
     &    num_multipolar),
     & map_placeholder5_to_deepbask(num_placeholder5_to_deepbask,
     &    num_deepbask),
     & map_placeholder5_to_deepng  (num_placeholder5_to_deepng  ,
     &    num_deepng),
     & map_placeholder5_to_deepLTS(num_placeholder5_to_deepLTS,
     &    num_deepLTS),
     &map_placeholder5_to_placeholder6(num_placeholder5_to_placeholder6,
     &    num_placeholder6),
     & map_placeholder5_to_L3pyr(num_placeholder5_to_L3pyr,num_L3pyr), 
     &map_placeholder6_to_placeholder5(num_placeholder6_to_placeholder5,
     &    num_placeholder5),
     &map_placeholder6_to_placeholder6(num_placeholder6_to_placeholder6,
     &    num_placeholder6),
     & map_L3pyr_to_L2pyr(num_L3pyr_to_L2pyr,
     &                             num_L2pyr), 
     & map_L3pyr_to_placeholder1(num_L3pyr_to_placeholder1,
     &                             num_placeholder1), 
     & map_L3pyr_to_placeholder2(num_L3pyr_to_placeholder2,
     &                             num_placeholder2),
     & map_L3pyr_to_placeholder3(num_L3pyr_to_placeholder3,
     &                             num_placeholder3),  
     & map_L3pyr_to_LOT(num_L3pyr_to_LOT,
     &                             num_LOT),
     & map_L3pyr_to_semilunar(num_L3pyr_to_semilunar,
     &                             num_semilunar),  
     & map_L3pyr_to_multipolar(num_L3pyr_to_multipolar,
     &                             num_multipolar),  
     & map_L3pyr_to_deepbask(num_L3pyr_to_deepbask,
     &                             num_deepbask), 
     & map_L3pyr_to_deepng  (num_L3pyr_to_deepng  ,
     &                             num_deepng  ), 
     & map_L3pyr_to_deepLTS(num_L3pyr_to_deepLTS,
     &                             num_deepLTS),
     & map_L3pyr_to_supVIP (num_L3pyr_to_supVIP ,
     &                             num_supVIP ),
     & map_L3pyr_to_placeholder5(num_L3pyr_to_placeholder5,
     &      num_placeholder5),
     & map_L3pyr_to_placeholder6(num_L3pyr_to_placeholder6,
     &      num_placeholder6)     
        INTEGER
     & map_L3pyr_to_L3pyr (num_L3pyr_to_L3pyr, num_L3pyr)

c Maps of synaptic compartments.  For simplicity, all contacts
c only made to one compartment.  Axoaxonic cells forced to contact 
c initial axon segments; other compartments will be listed in arrays.
! - except in original piriform model, no axoaxonic cells
        INTEGER 
     & com_L2pyr_to_L2pyr(num_L2pyr_to_L2pyr,
     &                           num_L2pyr),
     & com_L2pyr_to_placeholder1(num_L2pyr_to_placeholder1,  
     &                           num_placeholder1), 
     & com_L2pyr_to_supng  (num_L2pyr_to_supng  ,  
     &                           num_supng  ), 
     & com_L2pyr_to_placeholder2(num_L2pyr_to_placeholder2, 
     &                           num_placeholder2),
     & com_L2pyr_to_placeholder3(num_L2pyr_to_placeholder3,   
     &                           num_placeholder3),
     & com_L2pyr_to_LOT(num_L2pyr_to_LOT,
     &                           num_LOT),
     & com_L2pyr_to_semilunar(num_L2pyr_to_semilunar,
     &                           num_semilunar),  
     & com_L2pyr_to_multipolar(num_L2pyr_to_multipolar,
     &                           num_multipolar), 
     & com_L2pyr_to_deepbask(num_L2pyr_to_deepbask,
     &                           num_deepbask), 
     & com_L2pyr_to_deepLTS(num_L2pyr_to_deepLTS,
     &                           num_deepLTS), 
     & com_L2pyr_to_deepng  (num_L2pyr_to_deepng  ,
     &                             num_deepng  ), 
     & com_L2pyr_to_supVIP (num_L2pyr_to_supVIP ,
     &                           num_supVIP ), 
     & com_L2pyr_to_L3pyr(num_L2pyr_to_L3pyr,
     &                           num_L3pyr) 
              INTEGER
     & com_placeholder1_to_L2pyr(num_placeholder1_to_L2pyr,
     &                           num_L2pyr),  
     &com_placeholder1_to_placeholder1(num_placeholder1_to_placeholder1,
     &                           num_placeholder1), 
     & com_placeholder1_to_supng  (num_placeholder1_to_supng  ,
     &                           num_supng  ), 
     &com_placeholder1_to_placeholder2(num_placeholder1_to_placeholder2,
     &                           num_placeholder2),
     &com_placeholder1_to_placeholder3(num_placeholder1_to_placeholder3,
     &                           num_placeholder3),  
     & com_placeholder1_to_LOT(num_placeholder1_to_LOT,
     &                           num_LOT)  

          INTEGER
     & com_supng_to_L2pyr  (num_supng_to_L2pyr,
     &                         num_L2pyr),
     & com_supng_to_L3pyr (num_supng_to_L3pyr,
     &                         num_L3pyr),
     & com_supng_to_semilunar    (num_supng_to_semilunar   ,
     &                         num_semilunar   ),
     & com_supng_to_multipolar    (num_supng_to_multipolar   ,
     &                         num_multipolar   ),
     & com_supng_to_supng     (num_supng_to_supng    ,
     &                         num_supng    ),
     & com_supng_to_placeholder1   (num_supng_to_placeholder1  ,
     &                         num_placeholder1  ) 

          INTEGER
     & com_placeholder2_to_L2pyr(num_placeholder2_to_L2pyr,
     &                           num_L2pyr), 
     & com_placeholder2_to_LOT(num_placeholder2_to_LOT,
     &                           num_LOT)
           INTEGER
     & com_placeholder2_to_semilunar(num_placeholder2_to_semilunar,
     &                           num_semilunar),  
     &com_placeholder2_to_multipolar(num_placeholder2_to_multipolar,
     &                           num_multipolar), 
     & com_placeholder2_to_L3pyr(num_placeholder2_to_L3pyr,
     &                           num_L3pyr), 
     & com_placeholder3_to_L2pyr(num_placeholder3_to_L2pyr,
     &                           num_L2pyr),  
     &com_placeholder3_to_placeholder1(num_placeholder3_to_placeholder1,
     &                           num_placeholder1),  
     &com_placeholder3_to_placeholder2(num_placeholder3_to_placeholder2,
     &                           num_placeholder2), 
     &com_placeholder3_to_placeholder3(num_placeholder3_to_placeholder3,
     &                           num_placeholder3), 
     & com_placeholder3_to_LOT(num_placeholder3_to_LOT,
     &                           num_LOT), 
     & com_placeholder3_to_semilunar(num_placeholder3_to_semilunar,
     &                           num_semilunar),   
     &com_placeholder3_to_multipolar(num_placeholder3_to_multipolar,
     &                           num_multipolar),  
     & com_placeholder3_to_deepbask(num_placeholder3_to_deepbask,
     &                           num_deepbask), 
     & com_placeholder3_to_deepLTS(num_placeholder3_to_deepLTS,
     &                           num_deepLTS), 
     & com_placeholder3_to_supVIP (num_placeholder3_to_supVIP ,
     &                           num_supVIP ), 
     & com_placeholder3_to_L3pyr(num_placeholder3_to_L3pyr,
     &                           num_L3pyr), 
     & com_LOT_to_L2pyr(num_LOT_to_L2pyr,
     &                           num_L2pyr),
     & com_LOT_to_placeholder1(num_LOT_to_placeholder1,
     &                           num_placeholder1), 
     & com_LOT_to_placeholder2(num_LOT_to_placeholder2,
     &                           num_placeholder2)
                INTEGER
     & com_LOT_to_placeholder3(num_LOT_to_placeholder3,
     &                           num_placeholder3), 
     & com_LOT_to_LOT(num_LOT_to_LOT,
     &                           num_LOT),
     & com_LOT_to_semilunar(num_LOT_to_semilunar,
     &                           num_semilunar),  
     & com_LOT_to_multipolar(num_LOT_to_multipolar,
     &                           num_multipolar), 
     & com_LOT_to_deepbask(num_LOT_to_deepbask,
     &                           num_deepbask), 
     & com_LOT_to_deepng  (num_LOT_to_deepng  ,
     &                           num_deepng  ), 
     & com_LOT_to_deepLTS(num_LOT_to_deepLTS,
     &                           num_deepLTS),
     & com_LOT_to_supVIP (num_LOT_to_supVIP ,
     &                           num_supVIP ),
     & com_LOT_to_supng  (num_LOT_to_supng  ,
     &                           num_supng  ),
     & com_LOT_to_L3pyr(num_LOT_to_L3pyr,
     &                           num_L3pyr),
     & com_semilunar_to_L2pyr(num_semilunar_to_L2pyr,
     &                           num_L2pyr),   
     & com_semilunar_to_placeholder1(num_semilunar_to_placeholder1,
     &                           num_placeholder1),  
     & com_semilunar_to_placeholder2(num_semilunar_to_placeholder2,
     &                           num_placeholder2), 
     & com_semilunar_to_placeholder3(num_semilunar_to_placeholder3,
     &                           num_placeholder3), 
     & com_semilunar_to_LOT(num_semilunar_to_LOT,
     &                           num_LOT), 
     & com_semilunar_to_semilunar(num_semilunar_to_semilunar,
     &                           num_semilunar),   
     & com_semilunar_to_multipolar(num_semilunar_to_multipolar,
     &                           num_multipolar),  
     & com_semilunar_to_deepbask(num_semilunar_to_deepbask,
     &                           num_deepbask), 
     & com_semilunar_to_deepng  (num_semilunar_to_deepng  ,
     &                           num_deepng  ), 
     & com_semilunar_to_deepLTS(num_semilunar_to_deepLTS,
     &                           num_deepLTS),  
     & com_semilunar_to_supVIP (num_semilunar_to_supVIP ,
     &                           num_supVIP ),  
     & com_semilunar_to_L3pyr(num_semilunar_to_L3pyr,
     &                           num_L3pyr) 
              INTEGER
     & com_multipolar_to_L2pyr(num_multipolar_to_L2pyr,
     &                           num_L2pyr), 
     &com_multipolar_to_placeholder1(num_multipolar_to_placeholder1,
     &                           num_placeholder1),  
     &com_multipolar_to_placeholder2(num_multipolar_to_placeholder2,
     &                           num_placeholder2),   
     &com_multipolar_to_placeholder3(num_multipolar_to_placeholder3,
     &                           num_placeholder3),     
     & com_multipolar_to_LOT(num_multipolar_to_LOT,
     &                           num_LOT), 
     & com_multipolar_to_semilunar(num_multipolar_to_semilunar,
     &                           num_semilunar),   
     &com_multipolar_to_multipolar(num_multipolar_to_multipolar,
     &                           num_multipolar),     
     & com_multipolar_to_deepbask(num_multipolar_to_deepbask,
     &                           num_deepbask),  
     & com_multipolar_to_deepng  (num_multipolar_to_deepng  ,
     &                           num_deepng  ),  
     & com_multipolar_to_deepLTS(num_multipolar_to_deepLTS,
     &                           num_deepLTS),   
     & com_multipolar_to_supVIP (num_multipolar_to_supVIP ,
     &                           num_supVIP ),   
     & com_multipolar_to_L3pyr(num_multipolar_to_L3pyr,
     &                           num_L3pyr),  
     & com_deepbask_to_LOT(num_deepbask_to_LOT,
     &                           num_LOT), 
     & com_deepbask_to_semilunar(num_deepbask_to_semilunar,
     &                           num_semilunar),   
     & com_deepbask_to_multipolar(num_deepbask_to_multipolar,
     &                           num_multipolar),  
     & com_deepbask_to_deepbask(num_deepbask_to_deepbask,
     &                           num_deepbask), 
     & com_deepbask_to_deepng  (num_deepbask_to_deepng  ,
     &                           num_deepng  ), 
     & com_deepbask_to_deepLTS(num_deepbask_to_deepLTS,
     &                           num_deepLTS),  
     & com_deepbask_to_supVIP (num_deepbask_to_supVIP ,
     &                           num_supVIP ),  
     & com_deepbask_to_L2pyr(num_deepbask_to_L2pyr,
     &                           num_L2pyr), 
     & com_deepbask_to_L3pyr(num_deepbask_to_L3pyr,
     &                           num_L3pyr) 
            INTEGER
     & com_deepng_to_semilunar     (num_deepng_to_semilunar    ,
     &                           num_semilunar      ),
     & com_deepng_to_multipolar     (num_deepng_to_multipolar    ,
     &                           num_multipolar      ),
     & com_deepng_to_L2pyr  (num_deepng_to_L2pyr ,
     &                           num_L2pyr   ),
     & com_deepng_to_L3pyr  (num_deepng_to_L3pyr ,
     &                           num_L3pyr   ),
     & com_deepng_to_LOT  (num_deepng_to_LOT ,
     &                           num_LOT   ),
     & com_deepng_to_deepng     (num_deepng_to_deepng    ,
     &                           num_deepng      ),
     & com_deepng_to_deepbask   (num_deepng_to_deepbask  ,
     &                           num_deepbask    ) 
            INTEGER
     & com_deepLTS_to_L2pyr(num_deepLTS_to_L2pyr,
     &                           num_L2pyr), 
     & com_deepLTS_to_LOT(num_deepLTS_to_LOT,
     &                           num_LOT),
     & com_deepLTS_to_semilunar(num_deepLTS_to_semilunar,
     &                           num_semilunar), 
     & com_deepLTS_to_multipolar(num_deepLTS_to_multipolar,
     &                           num_multipolar),    
     & com_deepLTS_to_L3pyr(num_deepLTS_to_L3pyr,
     &                           num_L3pyr),

     & com_supVIP_to_L2pyr(num_supVIP_to_L2pyr,
     &                           num_L2pyr), 
     & com_supVIP_to_placeholder1(num_supVIP_to_placeholder1,
     &                           num_placeholder1),  
     & com_supVIP_to_placeholder2(num_supVIP_to_placeholder2,
     &                           num_placeholder2), 
     & com_supVIP_to_placeholder3(num_supVIP_to_placeholder3,
     &                           num_placeholder3), 
     & com_supVIP_to_supng (num_supVIP_to_supng ,
     &                           num_supng ), 
     & com_supVIP_to_LOT(num_supVIP_to_LOT,
     &                           num_LOT),
     & com_supVIP_to_semilunar(num_supVIP_to_semilunar,
     &                           num_semilunar),  
     & com_supVIP_to_multipolar(num_supVIP_to_multipolar,
     &                            num_multipolar), 
     & com_supVIP_to_deepbask(num_supVIP_to_deepbask,
     &                            num_deepbask), 
     & com_supVIP_to_deepLTS(num_supVIP_to_deepLTS,
     &                            num_deepLTS),  
     & com_supVIP_to_supVIP(num_supVIP_to_supVIP,
     &                            num_supVIP)  
           INTEGER
     & com_supVIP_to_L3pyr(num_supVIP_to_L3pyr,
     &                            num_L3pyr), 
     & com_placeholder5_to_L2pyr(num_placeholder5_to_L2pyr,
     &                            num_L2pyr),     
     & com_placeholder5_to_placeholder1(num_placeholder5_to_placeholder1
     &          ,                 num_placeholder1),    
     & com_placeholder5_to_supng  (num_placeholder5_to_supng  ,
     &                            num_supng  ),    
     & com_placeholder5_to_placeholder2(num_placeholder5_to_placeholder2
     &    ,num_placeholder2),
     & com_placeholder5_to_LOT(num_placeholder5_to_LOT,num_LOT),
     & com_placeholder5_to_semilunar(num_placeholder5_to_semilunar,
     &     num_semilunar),
     &com_placeholder5_to_multipolar(num_placeholder5_to_multipolar,
     &         num_multipolar),     
     & com_placeholder5_to_deepbask(num_placeholder5_to_deepbask,
     &     num_deepbask),
     &com_placeholder5_to_deepng(num_placeholder5_to_deepng,num_deepng), 
     & com_placeholder5_to_deepLTS(num_placeholder5_to_deepLTS,
     &     num_deepLTS),
     &com_placeholder5_to_placeholder6(num_placeholder5_to_placeholder6,
     &   num_placeholder6),
     & com_placeholder5_to_L3pyr(num_placeholder5_to_L3pyr,num_L3pyr), 
     &com_placeholder6_to_placeholder5(num_placeholder6_to_placeholder5,
     &   num_placeholder5),
     &com_placeholder6_to_placeholder6(num_placeholder6_to_placeholder6,
     &    num_placeholder6),
     & com_L3pyr_to_L2pyr(num_L3pyr_to_L2pyr,
     &                             num_L2pyr), 
     & com_L3pyr_to_placeholder1(num_L3pyr_to_placeholder1,
     &                             num_placeholder1), 
     & com_L3pyr_to_placeholder2(num_L3pyr_to_placeholder2,
     &                             num_placeholder2),
     & com_L3pyr_to_placeholder3(num_L3pyr_to_placeholder3,
     &                             num_placeholder3),  
     & com_L3pyr_to_LOT(num_L3pyr_to_LOT,
     &                             num_LOT),
     & com_L3pyr_to_semilunar(num_L3pyr_to_semilunar,
     &                             num_semilunar)  
              INTEGER
     & com_L3pyr_to_multipolar(num_L3pyr_to_multipolar,
     &                             num_multipolar),  
     & com_L3pyr_to_deepbask(num_L3pyr_to_deepbask,
     &                             num_deepbask), 
     & com_L3pyr_to_deepng  (num_L3pyr_to_deepng  ,
     &                             num_deepng  ), 
     & com_L3pyr_to_deepLTS(num_L3pyr_to_deepLTS,
     &                             num_deepLTS),
     & com_L3pyr_to_supVIP (num_L3pyr_to_supVIP ,
     &                             num_supVIP ),
     &com_L3pyr_to_placeholder5(num_L3pyr_to_placeholder5,
     &     num_placeholder5),
     &com_L3pyr_to_placeholder6(num_L3pyr_to_placeholder6,
     &     num_placeholder6),
     & com_L3pyr_to_L3pyr(num_L3pyr_to_L3pyr,
     &                             num_L3pyr)

        integer num_LOTstim_L2pyr, num_LOTstim_L3pyr,
     &          num_LOTstim_semilunar

c Entries in gjtable are cell a, compart. of cell a with gj,
c  cell b, compart. of cell b with gj; entries not repeated,
c which means that, for given cell being integrated, table
c must be searched through cols. 1 and 3.
       integer gjtable_L2pyr(totaxgj_L2pyr,4),
     &   gjtable_placeholder1  (totSDgj_placeholder1,4),
     &   gjtable_supng    (totSDgj_supng  ,4),
     &   gjtable_placeholder2  (1              ,4),
     &   gjtable_placeholder3   (totSDgj_placeholder3,4),
     &   gjtable_LOT(totaxgj_LOT,4),
     &   gjtable_semilunar   (totaxgj_semilunar,4),
     &   gjtable_multipolar   (totaxgj_multipolar,4),
     &   gjtable_L3pyr(totaxgj_L3pyr,4),
     &   gjtable_deepbask (totSDgj_deepbask,4),
     &   gjtable_deepng   (totSDgj_deepng  ,4),
     &   gjtable_deepLTS (1               ,4),
     &   gjtable_supVIP   (totSDgj_supVIP ,4),
     &   gjtable_placeholder5      (totaxgj_placeholder5,4),
     &   gjtable_placeholder6      (totaxgj_placeholder6,4) 

c define compartments on which gj can form
       INTEGER
     &table_axgjcompallow_L2pyr(num_axgjcompallow_L2pyr)
c    &          /74/,
     &          /73/, ! 28 Nov. 2005, move proximally, to get more inhib. control.
c Ectopics to L2pyr then go to #72, see
c   supergj.f
     &table_SDgjcompallow_placeholder1 (num_SDgjcompallow_placeholder1 )
     &          /3,4,16,17,29,30,42,43/,
     &table_SDgjcompallow_supng    (num_SDgjcompallow_supng    )
     &          /3,4,16,17,29,30,42,43/,
     &table_SDgjcompallow_placeholder3 (num_SDgjcompallow_placeholder3 )
     &          /3,4,16,17,29,30,42,43/,
     &table_axgjcompallow_LOT(num_axgjcompallow_LOT)
     &          /59/,
     &table_axgjcompallow_semilunar   (num_axgjcompallow_semilunar   )
c    &          /74/,
     &          /73/,
     &table_axgjcompallow_multipolar(num_axgjcompallow_multipolar  )
     &          /56/,
     &table_axgjcompallow_L3pyr(num_axgjcompallow_L3pyr)
c    &          /74/,
     &          /73/,
     &table_SDgjcompallow_deepbask (num_SDgjcompallow_deepbask )
     &          /3,4,16,17,29,30,42,43/,
     &table_SDgjcompallow_deepng   (num_SDgjcompallow_deepng   )
     &          /3,4,16,17,29,30,42,43/,
     &table_SDgjcompallow_supVIP   (num_SDgjcompallow_supVIP   )
     &          /3,4,16,17,29,30,42,43/,
     &table_axgjcompallow_placeholder5 (num_axgjcompallow_placeholder5)
     &          /137/,
c Ectopics to placeholder5 cells to #135
     &table_axgjcompallow_placeholder6 (num_axgjcompallow_placeholder6 )
     &          /57/


       real*8 field_sup, field_deep ! scalars to pass to subroutines
       real*8 field_sup_local(1), field_deep_local(1)  ! for mpi
       real*8 field_sup_global(numnodes), field_deep_global(numnodes) ! for mpi
       real*8 field_sup_tot, field_deep_tot  ! sums of global vectors

c Define tables used for computing dexp & GABA-B timecourse:
c dexptablesmall(i) = dexp(-z), i = int (z*1000.), 0<=z<=5.
c dexptablebig  (i) = dexp(-z), i = int (z*10.), 0<=z<=100.
        double precision:: dexptablesmall(0:5000)
        double precision::  dexptablebig  (0:1000)
        double precision:: otis_table (0:50000)
! if how_often = 50 and dt = .002, then otis_table structure
! corresponds to time steps of 0.1 ms, and it gives 5 s of data.

        real*8 noisepe_semilunar  ! noisepe_semilunar_save defined as parameter above
        real*8 gapcon_L2pyr
        real*8 z1ai, z1bi, z1ap, z1bp

c Define arrays, constants, for voltages, applied currents,
c synaptic conductances, random numbers, etc.

       double precision::
     &  V_L2pyr  (numcomp_L2pyr, num_L2pyr),
     &  V_placeholder1   (numcomp_placeholder1,  num_placeholder1),  
     &  V_supng     (numcomp_supng  ,  num_supng  ),  
     &  V_placeholder2   (numcomp_placeholder2,  num_placeholder2), 
     &  V_placeholder3    (numcomp_placeholder3,   num_placeholder3), 
     &  V_LOT (numcomp_LOT,num_LOT),
     &  V_semilunar    (numcomp_semilunar,   num_semilunar),  
     &  V_multipolar    (numcomp_multipolar,   num_multipolar), 
     &  V_L3pyr (numcomp_L3pyr,num_L3pyr),
     &  V_deepbask  (numcomp_deepbask, num_deepbask),
     &  V_deepng    (numcomp_deepng  , num_deepng  ),
     &  V_deepLTS  (numcomp_deepLTS, num_deepLTS),
     &  V_supVIP    (numcomp_supVIP ,  num_supVIP ),
     &  V_placeholder5  (numcomp_placeholder5,      num_placeholder5),   
     &  V_placeholder6  (numcomp_placeholder6,      num_placeholder6) 

       double precision::
     &  curr_L2pyr   (numcomp_L2pyr, num_L2pyr),
     &  curr_placeholder1   (numcomp_placeholder1,  num_placeholder1),  
     &  curr_supng      (numcomp_supng  ,  num_supng  ),  
     &  curr_placeholder2    (numcomp_placeholder2,  num_placeholder2), 
     &  curr_placeholder3    (numcomp_placeholder3,   num_placeholder3), 
     &  curr_LOT  (numcomp_LOT,num_LOT),
     &  curr_semilunar     (numcomp_semilunar,   num_semilunar),  
     &  curr_multipolar   (numcomp_multipolar,   num_multipolar), 
     &  curr_L3pyr  (numcomp_L3pyr,num_L3pyr),
     &  curr_deepbask   (numcomp_deepbask, num_deepbask),
     &  curr_deepng     (numcomp_deepng  , num_deepng  ),
     &  curr_deepLTS   (numcomp_deepLTS, num_deepLTS),
     &  curr_supVIP     (numcomp_supVIP ,  num_supVIP ),
     &  curr_placeholder5 (numcomp_placeholder5,      num_placeholder5),   
     &  curr_placeholder6  (numcomp_placeholder6,      num_placeholder6) 

       double precision::
     & gAMPA_L2pyr   (numcomp_L2pyr, num_L2pyr),
     & gAMPA_placeholder1    (numcomp_placeholder1,  num_placeholder1),  
     & gAMPA_supng      (numcomp_supng  ,  num_supng  ),  
     & gAMPA_placeholder2    (numcomp_placeholder2,  num_placeholder2), 
     & gAMPA_placeholder3    (numcomp_placeholder3,   num_placeholder3), 
     & gAMPA_LOT  (numcomp_LOT,num_LOT),
     & gAMPA_semilunar     (numcomp_semilunar,   num_semilunar),  
     & gAMPA_multipolar    (numcomp_multipolar,   num_multipolar), 
     & gAMPA_L3pyr  (numcomp_L3pyr,num_L3pyr),
     & gAMPA_deepbask   (numcomp_deepbask, num_deepbask),
     & gAMPA_deepng     (numcomp_deepng  , num_deepng  ),
     & gAMPA_deepLTS   (numcomp_deepLTS, num_deepLTS),
     & gAMPA_supVIP     (numcomp_supVIP ,  num_supVIP ),
     & gAMPA_placeholder5 (numcomp_placeholder5,      num_placeholder5),   
     & gAMPA_placeholder6  (numcomp_placeholder6,      num_placeholder6) 

       double precision::
     & gNMDA_L2pyr   (numcomp_L2pyr, num_L2pyr),
     & gNMDA_placeholder1    (numcomp_placeholder1,  num_placeholder1),  
     & gNMDA_supng      (numcomp_supng  ,  num_supng  ),  
     & gNMDA_placeholder2    (numcomp_placeholder2,  num_placeholder2), 
     & gNMDA_placeholder3    (numcomp_placeholder3,   num_placeholder3), 
     & gNMDA_LOT  (numcomp_LOT,num_LOT),
     & gNMDA_semilunar     (numcomp_semilunar,   num_semilunar),  
     & gNMDA_multipolar     (numcomp_multipolar,  num_multipolar), 
     & gNMDA_L3pyr  (numcomp_L3pyr,num_L3pyr),
     & gNMDA_deepbask   (numcomp_deepbask, num_deepbask),
     & gNMDA_deepng     (numcomp_deepng  , num_deepng  ),
     & gNMDA_deepLTS   (numcomp_deepLTS, num_deepLTS),
     & gNMDA_supVIP     (numcomp_supVIP ,  num_supVIP ),
     & gNMDA_placeholder5 (numcomp_placeholder5,      num_placeholder5),   
     & gNMDA_placeholder6 (numcomp_placeholder6,      num_placeholder6) 

       double precision::
     & gGABA_A_L2pyr (numcomp_L2pyr, num_L2pyr),
     & gGABA_A_placeholder1  (numcomp_placeholder1,  num_placeholder1),  
     & gGABA_A_supng    (numcomp_supng  ,  num_supng  ),  
     & gGABA_A_placeholder2  (numcomp_placeholder2,  num_placeholder2), 
     & gGABA_A_placeholder3  (numcomp_placeholder3,   num_placeholder3), 
     & gGABA_A_LOT(numcomp_LOT,num_LOT),
     & gGABA_A_semilunar   (numcomp_semilunar,   num_semilunar),  
     & gGABA_A_multipolar  (numcomp_multipolar,   num_multipolar), 
     & gGABA_A_L3pyr(numcomp_L3pyr,num_L3pyr),
     & gGABA_A_deepbask (numcomp_deepbask, num_deepbask),
     & gGABA_A_deepng   (numcomp_deepng  , num_deepng  ),
     & gGABA_A_deepLTS (numcomp_deepLTS, num_deepLTS),
     & gGABA_A_supVIP   (numcomp_supVIP ,  num_supVIP ),
     & gGABA_A_placeholder5(numcomp_placeholder5,    num_placeholder5),   
     & gGABA_A_placeholder6(numcomp_placeholder6,    num_placeholder6) 

       double precision::
     & gGABA_B_L2pyr (numcomp_L2pyr, num_L2pyr),
     & gGABA_B_LOT(numcomp_LOT,num_LOT),
     & gGABA_B_semilunar  (numcomp_semilunar,   num_semilunar),  
c    & gGABA_B_multipolar (numcomp_multipolar,   num_multipolar), 
c Perhaps modify integrate_multipolar later to deal with GABA-B
     & gGABA_B_L3pyr(numcomp_L3pyr,num_L3pyr),
     & gGABA_B_placeholder5 (numcomp_placeholder5,    num_placeholder5),   
     & gGABA_B_placeholder6 (numcomp_placeholder6,    num_placeholder6) 

! define membrane and Ca state variables that must be passed
! to subroutines
       real*8  chi_L2pyr(numcomp_L2pyr,num_L2pyr)
       real*8  mnaf_L2pyr(numcomp_L2pyr,num_L2pyr),
     & mnap_L2pyr(numcomp_L2pyr,num_L2pyr),
     x hnaf_L2pyr(numcomp_L2pyr,num_L2pyr),
     x mkdr_L2pyr(numcomp_L2pyr,num_L2pyr),
     x mka_L2pyr(numcomp_L2pyr,num_L2pyr),
     x hka_L2pyr(numcomp_L2pyr,num_L2pyr),
     x mk2_L2pyr(numcomp_L2pyr,num_L2pyr), 
     x hk2_L2pyr(numcomp_L2pyr,num_L2pyr),
     x mkm_L2pyr(numcomp_L2pyr,num_L2pyr),
     x mkc_L2pyr(numcomp_L2pyr,num_L2pyr),
     x mkahp_L2pyr(numcomp_L2pyr,num_L2pyr),
     x mcat_L2pyr(numcomp_L2pyr,num_L2pyr),
     x hcat_L2pyr(numcomp_L2pyr,num_L2pyr),
     x mcal_L2pyr(numcomp_L2pyr,num_L2pyr),
     x mar_L2pyr(numcomp_L2pyr,num_L2pyr)

       real*8  chi_placeholder1 (numcomp_placeholder1 ,num_placeholder1)
       real*8  mnaf_placeholder1(numcomp_placeholder1,num_placeholder1),
     & mnap_placeholder1 (numcomp_placeholder1 ,num_placeholder1 ),
     x hnaf_placeholder1 (numcomp_placeholder1 ,num_placeholder1 ),
     x mkdr_placeholder1 (numcomp_placeholder1 ,num_placeholder1 ),
     x mka_placeholder1 (numcomp_placeholder1 ,num_placeholder1 ),
     x hka_placeholder1 (numcomp_placeholder1 ,num_placeholder1 ),
     x mk2_placeholder1 (numcomp_placeholder1 ,num_placeholder1 ), 
     x hk2_placeholder1 (numcomp_placeholder1 ,num_placeholder1 ),
     x mkm_placeholder1 (numcomp_placeholder1 ,num_placeholder1 ),
     x mkc_placeholder1 (numcomp_placeholder1 ,num_placeholder1 ),
     x mkahp_placeholder1 (numcomp_placeholder1 ,num_placeholder1 ),
     x mcat_placeholder1 (numcomp_placeholder1 ,num_placeholder1 ),
     x hcat_placeholder1 (numcomp_placeholder1 ,num_placeholder1 ),
     x mcal_placeholder1 (numcomp_placeholder1 ,num_placeholder1 ),
     x mar_placeholder1 (numcomp_placeholder1 ,num_placeholder1 )

       real*8  chi_supng (numcomp_supng ,num_supng )
       real*8  mnaf_supng (numcomp_supng ,num_supng ),
     & mnap_supng (numcomp_supng ,num_supng ),
     x hnaf_supng (numcomp_supng ,num_supng ),
     x mkdr_supng (numcomp_supng ,num_supng ),
     x mka_supng (numcomp_supng ,num_supng ),
     x hka_supng (numcomp_supng ,num_supng ),
     x mk2_supng (numcomp_supng ,num_supng ), 
     x hk2_supng (numcomp_supng ,num_supng ),
     x mkm_supng (numcomp_supng ,num_supng ),
     x mkc_supng (numcomp_supng ,num_supng ),
     x mkahp_supng (numcomp_supng ,num_supng ),
     x mcat_supng (numcomp_supng ,num_supng ),
     x hcat_supng (numcomp_supng ,num_supng ),
     x mcal_supng (numcomp_supng ,num_supng ),
     x mar_supng (numcomp_supng ,num_supng )

       real*8  chi_placeholder2 (numcomp_placeholder2 ,num_placeholder2)
       real*8  mnaf_placeholder2(numcomp_placeholder2,num_placeholder2),
     & mnap_placeholder2 (numcomp_placeholder2 ,num_placeholder2 ),
     x hnaf_placeholder2 (numcomp_placeholder2 ,num_placeholder2 ),
     x mkdr_placeholder2 (numcomp_placeholder2 ,num_placeholder2 ),
     x mka_placeholder2 (numcomp_placeholder2 ,num_placeholder2 ),
     x hka_placeholder2 (numcomp_placeholder2 ,num_placeholder2 ),
     x mk2_placeholder2 (numcomp_placeholder2 ,num_placeholder2 ), 
     x hk2_placeholder2 (numcomp_placeholder2 ,num_placeholder2 ),
     x mkm_placeholder2 (numcomp_placeholder2 ,num_placeholder2 ),
     x mkc_placeholder2 (numcomp_placeholder2 ,num_placeholder2 ),
     x mkahp_placeholder2 (numcomp_placeholder2 ,num_placeholder2 ),
     x mcat_placeholder2 (numcomp_placeholder2 ,num_placeholder2 ),
     x hcat_placeholder2 (numcomp_placeholder2 ,num_placeholder2 ),
     x mcal_placeholder2 (numcomp_placeholder2 ,num_placeholder2 ),
     x mar_placeholder2 (numcomp_placeholder2 ,num_placeholder2 )

       real*8  chi_placeholder3(numcomp_placeholder3,num_placeholder3)
       real*8  mnaf_placeholder3(numcomp_placeholder3,num_placeholder3),
     & mnap_placeholder3(numcomp_placeholder3,num_placeholder3),
     x hnaf_placeholder3(numcomp_placeholder3,num_placeholder3),
     x mkdr_placeholder3(numcomp_placeholder3,num_placeholder3),
     x mka_placeholder3(numcomp_placeholder3,num_placeholder3),
     x hka_placeholder3(numcomp_placeholder3,num_placeholder3),
     x mk2_placeholder3(numcomp_placeholder3,num_placeholder3), 
     x hk2_placeholder3(numcomp_placeholder3,num_placeholder3),
     x mkm_placeholder3(numcomp_placeholder3,num_placeholder3),
     x mkc_placeholder3(numcomp_placeholder3,num_placeholder3),
     x mkahp_placeholder3(numcomp_placeholder3,num_placeholder3),
     x mcat_placeholder3(numcomp_placeholder3,num_placeholder3),
     x hcat_placeholder3(numcomp_placeholder3,num_placeholder3),
     x mcal_placeholder3(numcomp_placeholder3,num_placeholder3),
     x mar_placeholder3(numcomp_placeholder3,num_placeholder3)

      real*8  chi_LOT(numcomp_LOT,num_LOT)
      real*8  mnaf_LOT(numcomp_LOT,num_LOT),
     & mnap_LOT(numcomp_LOT,num_LOT),
     x hnaf_LOT(numcomp_LOT,num_LOT),
     x mkdr_LOT(numcomp_LOT,num_LOT),
     x mka_LOT(numcomp_LOT,num_LOT),
     x hka_LOT(numcomp_LOT,num_LOT),
     x mk2_LOT(numcomp_LOT,num_LOT), 
     x hk2_LOT(numcomp_LOT,num_LOT),
     x mkm_LOT(numcomp_LOT,num_LOT),
     x mkc_LOT(numcomp_LOT,num_LOT),
     x mkahp_LOT(numcomp_LOT,num_LOT),
     x mcat_LOT(numcomp_LOT,num_LOT),
     x hcat_LOT(numcomp_LOT,num_LOT),
     x mcal_LOT(numcomp_LOT,num_LOT),
     x mar_LOT(numcomp_LOT,num_LOT)


       real*8  chi_semilunar(numcomp_semilunar,num_semilunar)
       real*8  mnaf_semilunar(numcomp_semilunar,num_semilunar),
     & mnap_semilunar(numcomp_semilunar,num_semilunar),
     x hnaf_semilunar(numcomp_semilunar,num_semilunar),
     x mkdr_semilunar(numcomp_semilunar,num_semilunar),
     x mka_semilunar(numcomp_semilunar,num_semilunar),
     x hka_semilunar(numcomp_semilunar,num_semilunar),
     x mk2_semilunar(numcomp_semilunar,num_semilunar), 
     x hk2_semilunar(numcomp_semilunar,num_semilunar),
     x mkm_semilunar(numcomp_semilunar,num_semilunar),
     x mkc_semilunar(numcomp_semilunar,num_semilunar),
     x mkahp_semilunar(numcomp_semilunar,num_semilunar),
     x mcat_semilunar(numcomp_semilunar,num_semilunar),
     x hcat_semilunar(numcomp_semilunar,num_semilunar),
     x mcal_semilunar(numcomp_semilunar,num_semilunar),
     x mar_semilunar(numcomp_semilunar,num_semilunar)

       real*8  chi_multipolar(numcomp_multipolar,num_multipolar)
       real*8  mnaf_multipolar(numcomp_multipolar,num_multipolar),
     & mnap_multipolar(numcomp_multipolar,num_multipolar),
     x hnaf_multipolar(numcomp_multipolar,num_multipolar),
     x mkdr_multipolar(numcomp_multipolar,num_multipolar),
     x mka_multipolar(numcomp_multipolar,num_multipolar),
     x hka_multipolar(numcomp_multipolar,num_multipolar),
     x mk2_multipolar(numcomp_multipolar,num_multipolar), 
     x hk2_multipolar(numcomp_multipolar,num_multipolar),
     x mkm_multipolar(numcomp_multipolar,num_multipolar),
     x mkc_multipolar(numcomp_multipolar,num_multipolar),
     x mkahp_multipolar(numcomp_multipolar,num_multipolar),
     x mcat_multipolar(numcomp_multipolar,num_multipolar),
     x hcat_multipolar(numcomp_multipolar,num_multipolar),
     x mcal_multipolar(numcomp_multipolar,num_multipolar),
     x mar_multipolar(numcomp_multipolar,num_multipolar)

       real*8  chi_L3pyr(numcomp_L3pyr,num_L3pyr)
       real*8  mnaf_L3pyr(numcomp_L3pyr,num_L3pyr),
     & mnap_L3pyr(numcomp_L3pyr,num_L3pyr),
     x hnaf_L3pyr(numcomp_L3pyr,num_L3pyr),
     x mkdr_L3pyr(numcomp_L3pyr,num_L3pyr),
     x mka_L3pyr(numcomp_L3pyr,num_L3pyr),
     x hka_L3pyr(numcomp_L3pyr,num_L3pyr),
     x mk2_L3pyr(numcomp_L3pyr,num_L3pyr), 
     x hk2_L3pyr(numcomp_L3pyr,num_L3pyr),
     x mkm_L3pyr(numcomp_L3pyr,num_L3pyr),
     x mkc_L3pyr(numcomp_L3pyr,num_L3pyr),
     x mkahp_L3pyr(numcomp_L3pyr,num_L3pyr),
     x mcat_L3pyr(numcomp_L3pyr,num_L3pyr),
     x hcat_L3pyr(numcomp_L3pyr,num_L3pyr),
     x mcal_L3pyr(numcomp_L3pyr,num_L3pyr),
     x mar_L3pyr(numcomp_L3pyr,num_L3pyr)

       real*8  chi_deepbask(numcomp_deepbask,num_deepbask)
       real*8  mnaf_deepbask(numcomp_deepbask,num_deepbask),
     & mnap_deepbask(numcomp_deepbask,num_deepbask),
     x hnaf_deepbask(numcomp_deepbask,num_deepbask),
     x mkdr_deepbask(numcomp_deepbask,num_deepbask),
     x mka_deepbask(numcomp_deepbask,num_deepbask),
     x hka_deepbask(numcomp_deepbask,num_deepbask),
     x mk2_deepbask(numcomp_deepbask,num_deepbask), 
     x hk2_deepbask(numcomp_deepbask,num_deepbask),
     x mkm_deepbask(numcomp_deepbask,num_deepbask),
     x mkc_deepbask(numcomp_deepbask,num_deepbask),
     x mkahp_deepbask(numcomp_deepbask,num_deepbask),
     x mcat_deepbask(numcomp_deepbask,num_deepbask),
     x hcat_deepbask(numcomp_deepbask,num_deepbask),
     x mcal_deepbask(numcomp_deepbask,num_deepbask),
     x mar_deepbask(numcomp_deepbask,num_deepbask)

       real*8  chi_deepng(numcomp_deepng,num_deepng)
       real*8  mnaf_deepng(numcomp_deepng,num_deepng),
     & mnap_deepng(numcomp_deepng,num_deepng),
     x hnaf_deepng(numcomp_deepng,num_deepng),
     x mkdr_deepng(numcomp_deepng,num_deepng),
     x mka_deepng(numcomp_deepng,num_deepng),
     x hka_deepng(numcomp_deepng,num_deepng),
     x mk2_deepng(numcomp_deepng,num_deepng), 
     x hk2_deepng(numcomp_deepng,num_deepng),
     x mkm_deepng(numcomp_deepng,num_deepng),
     x mkc_deepng(numcomp_deepng,num_deepng),
     x mkahp_deepng(numcomp_deepng,num_deepng),
     x mcat_deepng(numcomp_deepng,num_deepng),
     x hcat_deepng(numcomp_deepng,num_deepng),
     x mcal_deepng(numcomp_deepng,num_deepng),
     x mar_deepng(numcomp_deepng,num_deepng)

       real*8  chi_deepLTS(numcomp_deepLTS,num_deepLTS)
       real*8  mnaf_deepLTS(numcomp_deepLTS,num_deepLTS),
     & mnap_deepLTS(numcomp_deepLTS,num_deepLTS),
     x hnaf_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mkdr_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mka_deepLTS(numcomp_deepLTS,num_deepLTS),
     x hka_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mk2_deepLTS(numcomp_deepLTS,num_deepLTS), 
     x hk2_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mkm_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mkc_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mkahp_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mcat_deepLTS(numcomp_deepLTS,num_deepLTS),
     x hcat_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mcal_deepLTS(numcomp_deepLTS,num_deepLTS),
     x mar_deepLTS(numcomp_deepLTS,num_deepLTS)

       real*8  chi_supVIP(numcomp_supVIP,num_supVIP)
       real*8  mnaf_supVIP(numcomp_supVIP,num_supVIP),
     & mnap_supVIP(numcomp_supVIP,num_supVIP),
     x hnaf_supVIP(numcomp_supVIP,num_supVIP),
     x mkdr_supVIP(numcomp_supVIP,num_supVIP),
     x mka_supVIP(numcomp_supVIP,num_supVIP),
     x hka_supVIP(numcomp_supVIP,num_supVIP),
     x mk2_supVIP(numcomp_supVIP,num_supVIP), 
     x hk2_supVIP(numcomp_supVIP,num_supVIP),
     x mkm_supVIP(numcomp_supVIP,num_supVIP),
     x mkc_supVIP(numcomp_supVIP,num_supVIP),
     x mkahp_supVIP(numcomp_supVIP,num_supVIP),
     x mcat_supVIP(numcomp_supVIP,num_supVIP),
     x hcat_supVIP(numcomp_supVIP,num_supVIP),
     x mcal_supVIP(numcomp_supVIP,num_supVIP),
     x mar_supVIP(numcomp_supVIP,num_supVIP)

       real*8  chi_placeholder5(numcomp_placeholder5,num_placeholder5)
       real*8  mnaf_placeholder5(numcomp_placeholder5,num_placeholder5),
     & mnap_placeholder5(numcomp_placeholder5,num_placeholder5),
     x hnaf_placeholder5(numcomp_placeholder5,num_placeholder5),
     x mkdr_placeholder5(numcomp_placeholder5,num_placeholder5),
     x mka_placeholder5(numcomp_placeholder5,num_placeholder5),
     x hka_placeholder5(numcomp_placeholder5,num_placeholder5),
     x mk2_placeholder5(numcomp_placeholder5,num_placeholder5), 
     x hk2_placeholder5(numcomp_placeholder5,num_placeholder5),
     x mkm_placeholder5(numcomp_placeholder5,num_placeholder5),
     x mkc_placeholder5(numcomp_placeholder5,num_placeholder5),
     x mkahp_placeholder5(numcomp_placeholder5,num_placeholder5),
     x mcat_placeholder5(numcomp_placeholder5,num_placeholder5),
     x hcat_placeholder5(numcomp_placeholder5,num_placeholder5),
     x mcal_placeholder5(numcomp_placeholder5,num_placeholder5),
     x mar_placeholder5(numcomp_placeholder5,num_placeholder5)

       real*8  chi_placeholder6(numcomp_placeholder6,num_placeholder6)
       real*8  mnaf_placeholder6(numcomp_placeholder6,num_placeholder6),
     & mnap_placeholder6(numcomp_placeholder6,num_placeholder6),
     x hnaf_placeholder6(numcomp_placeholder6,num_placeholder6),
     x mkdr_placeholder6(numcomp_placeholder6,num_placeholder6),
     x mka_placeholder6(numcomp_placeholder6,num_placeholder6),
     x hka_placeholder6(numcomp_placeholder6,num_placeholder6),
     x mk2_placeholder6(numcomp_placeholder6,num_placeholder6), 
     x hk2_placeholder6(numcomp_placeholder6,num_placeholder6),
     x mkm_placeholder6(numcomp_placeholder6,num_placeholder6),
     x mkc_placeholder6(numcomp_placeholder6,num_placeholder6),
     x mkahp_placeholder6(numcomp_placeholder6,num_placeholder6),
     x mcat_placeholder6(numcomp_placeholder6,num_placeholder6),
     x hcat_placeholder6(numcomp_placeholder6,num_placeholder6),
     x mcal_placeholder6(numcomp_placeholder6,num_placeholder6),
     x mar_placeholder6(numcomp_placeholder6,num_placeholder6)

       double precision
     &    ranvec_L2pyr  (num_L2pyr),
     &    ranvec_placeholder1   (num_placeholder1),  
     &    ranvec_supng     (num_supng  ),  
     &    ranvec_placeholder2   (num_placeholder2), 
     &    ranvec_placeholder3    (num_placeholder3), 
     &    ranvec_LOT (num_LOT),
     &    ranvec_semilunar    (num_semilunar),  
     &    ranvec_multipolar      (num_multipolar  ), 
     &    ranvec_L3pyr (num_L3pyr),
     &    ranvec_deepbask  (num_deepbask),
     &    ranvec_deepng    (num_deepng  ),
     &    ranvec_deepLTS  (num_deepLTS),
     &    ranvec_supVIP    (num_supVIP ),
     &    ranvec_placeholder5       (num_placeholder5),   
     &    ranvec_placeholder6       (num_placeholder6),
     &    seed /137.d0/

c Define arrays for distal axon voltages which will be shared
c between nodes, and for axonal sites of possible gj
         double precision::
     &  distal_axon_L2pyr  (maxcellspernode),
     &  ldistal_axon_L2pyr (num_L2pyr), ! use for outtime
     &  distal_axon_supintern (maxcellspernode),
     &  ldistal_axon_placeholder1   (num_placeholder1  ),
     &  ldistal_axon_placeholder2   (num_placeholder2  ),
     &  ldistal_axon_placeholder3    (num_placeholder3   ),
     &  ldistal_axon_supng     (num_supng    ),
     &  ldistal_axon_supVIP    (num_supVIP   )

         double precision::
     &  distal_axon_LOT (maxcellspernode),
     &  ldistal_axon_LOT(num_LOT),
     &  distal_axon_semilunar    (maxcellspernode),
     &  ldistal_axon_semilunar   (num_semilunar),
     &  distal_axon_multipolar      (maxcellspernode),
     &  ldistal_axon_multipolar     (num_multipolar  ),
     &  distal_axon_L3pyr (maxcellspernode),
     &  ldistal_axon_L3pyr(num_L3pyr),
     &  distal_axon_deepintern(maxcellspernode),
     &  ldistal_axon_deepbask (num_deepbask  ),
     &  ldistal_axon_deepLTS (num_deepLTS  ),
     &  ldistal_axon_deepng   (num_deepng    ),
     &  distal_axon_placeholder5       (maxcellspernode),
     &  ldistal_axon_placeholder5      (num_placeholder5),
     &  distal_axon_placeholder6       (maxcellspernode),
     &  ldistal_axon_placeholder6      (num_placeholder6),
!    Communication will be complicated, however, because - say - a semilunar
!   will have to communicate only the semilunar axons it has integrated.
     &  distal_axon_global    (numnodes  * maxcellspernode)
! distal_axon_global will be concatenation of individual
! distal_axon vectors       


         double precision::
     &  outtime_L2pyr  (5000, num_L2pyr),
     &  outtime_placeholder1   (5000, num_placeholder1), 
     &  outtime_supng     (5000, num_supng  ), 
     &  outtime_placeholder2   (5000, num_placeholder2), 
     &  outtime_placeholder3    (5000, num_placeholder3),   
     &  outtime_LOT (5000, num_LOT), 
     &  outtime_semilunar    (5000, num_semilunar), 
     &  outtime_multipolar      (5000, num_multipolar  ),  
     &  outtime_L3pyr (5000, num_L3pyr),
     &  outtime_deepbask  (5000, num_deepbask),
     &  outtime_deepng    (5000, num_deepng  ),
     &  outtime_deepLTS  (5000, num_deepLTS),
     &  outtime_supVIP    (5000, num_supVIP ), 
     &  outtime_placeholder5       (5000, num_placeholder5),      
     &  outtime_placeholder6       (5000, num_placeholder6)       

         INTEGER
     &  outctr_L2pyr  (num_L2pyr), 
     &  outctr_placeholder1   (num_placeholder1), 
     &  outctr_supng     (num_supng  ), 
     &  outctr_placeholder2   (num_placeholder2),
     &  outctr_placeholder3    (num_placeholder3),
     &  outctr_LOT (num_LOT),
     &  outctr_semilunar    (num_semilunar), 
     &  outctr_multipolar      (num_multipolar  ),
     &  outctr_L3pyr (num_L3pyr),
     &  outctr_deepbask  (num_deepbask),
     &  outctr_deepng    (num_deepng  ),
     &  outctr_deepLTS  (num_deepLTS),
     &  outctr_supVIP    (num_supVIP ),
     &  outctr_placeholder5       (num_placeholder5), 
     &  outctr_placeholder6       (num_placeholder6)

        CHARACTER(LEN=12) nodecell(0:numnodes-1) ! will define which cell type is to be handled by each node

        INTEGER place(0:numnodes-1)  ! this will define whether a node is 1st, 2nd... in the set of nodes
! used by a given type of cell

        integer initialize, firstcell, lastcell ! used in integration calls 
        integer ictr, ioffset

       REAL*8 gettime, time1, time2, time, timtot
       REAL*8 presyntime, delta, dexparg, dexparg1, dexparg2
       INTEGER thisno, display /0/, O
       REAL*8 z, z1, z2, outrcd(20), z3, z4, z3a, z4a, z5, z6, z7
       REAL*8 z10, z11, z12, z13, z14, z10a, z10b
       REAL*8 zz
       INTEGER i, j, k, L, k0, m

       double precision rel_axonshift_semilunar /0.d0/
       double precision rel_axonshift_L2pyr /0.d0/
       double precision rel_axonshift_L3pyr /0.d0/

c START EXECUTION PHASE
          include 'mpif.h'
          call mpi_init (info)
          call mpi_comm_rank(mpi_comm_world, thisno, info)
          call mpi_comm_size(mpi_comm_world, nodes , info)
          time1 = gettime()

c intialize outctr arrays
           do i = 1, num_L2pyr
        outctr_L2pyr  (i) = 0
           end do
           do i = 1, num_placeholder1  
        outctr_placeholder1 (i) = 0
           end do
           do i = 1, num_supng
        outctr_supng     (i) = 0
           end do
           do i = 1, num_placeholder2
        outctr_placeholder2 (i) = 0
           end do
           do i = 1, num_placeholder3
        outctr_placeholder3  (i) = 0
           end do
           do i = 1, num_LOT
        outctr_LOT (i) = 0
           end do
           do i = 1, num_semilunar
        outctr_semilunar  (i) = 0
           end do
           do i = 1, num_multipolar  
        outctr_multipolar    (i) = 0
           end do
           do i = 1, num_L3pyr
        outctr_L3pyr (i) = 0
           end do
           do i = 1, num_deepbask
        outctr_deepbask  (i) = 0
           end do
           do i = 1, num_deepng
        outctr_deepng    (i) = 0
           end do
           do i = 1, num_deepLTS
        outctr_deepLTS  (i) = 0
           end do
           do i = 1, num_supVIP
        outctr_supVIP    (i) = 0
           end do
           do i = 1, num_placeholder5
        outctr_placeholder5  (i) = 0
           end do
           do i = 1, num_placeholder6
        outctr_placeholder6  (i) = 0
           end do



c Define which cell type is handled by each processor
           nodecell(0) = 'L2pyr       '
           nodecell(1) = 'L2pyr       '
           nodecell(2) = 'supintern   '
           nodecell(3) = 'LOT         '
           nodecell(4) = 'semilunar   '
           nodecell(5) = 'multipolar  '
           nodecell(6) = 'L3pyr       '
           nodecell(7) = 'deepintern  '
           nodecell(8) = 'placeholder5'
           nodecell(9) = 'placeholder6'
          if (thisno.eq.0) then
            do i = 0, numnodes - 1
              write(6,786) i, nodecell(i)
786           format(i5,a12)
            end do
          end if

c Define "rank" of nodes assigned to each cell-type - will
c be used in figuring out how to partition the cells.
           place( 0) = 1  ! L2pyr: 1
           place( 1) = 2  ! L2pyr: 2
           place( 2) = 1  ! supintern  
           place( 3) = 1  ! LOT   
           place( 4) = 1  ! semilunar    
           place( 5) = 1  ! multipolar      
           place( 6) = 1  ! L3pyr 
           place( 7) = 1  ! deepintern
           place( 8) = 1  ! placeholder5       
           place( 9) = 1  ! placeholder6       

         do i = 1, 5000
           do j = 1, num_L2pyr
        outtime_L2pyr(i,j)             = -1.d5
           end do ! j
           do j = 1, num_placeholder1  
        outtime_placeholder1(i,j)              = -1.d5
           end do ! j
           do j = 1, num_supng    
        outtime_supng  (i,j)              = -1.d5
           end do ! j
           do j = 1, num_placeholder2  
        outtime_placeholder2(i,j)              = -1.d5
           end do ! j
           do j = 1, num_placeholder3   
        outtime_placeholder3(i,j)               = -1.d5
           end do ! j
           do j = 1, num_LOT
        outtime_LOT(i,j)            = -1.d5 
           end do ! j
           do j = 1, num_semilunar   
        outtime_semilunar(i,j)               = -1.d5
           end do ! j
           do j = 1, num_multipolar     
        outtime_multipolar(i,j)               = -1.d5
           end do ! j
           do j = 1, num_L3pyr   
        outtime_L3pyr(i,j)            = -1.d5
           end do ! j
           do j = 1, num_deepbask    
        outtime_deepbask(i,j)             = -1.d5
           end do ! j
           do j = 1, num_deepng      
        outtime_deepng  (i,j)             = -1.d5
           end do ! j
           do j = 1, num_deepLTS    
        outtime_deepLTS(i,j)             = -1.d5
           end do ! j
           do j = 1, num_supVIP      
        outtime_supVIP (i,j)              = -1.d5
           end do ! j
           do j = 1, num_placeholder5         
        outtime_placeholder5(i,j)                = -1.d5
           end do ! j
           do j = 1, num_placeholder6         
        outtime_placeholder6(i,j)                  = -1.d5
           end do ! j
         end do ! do i

c         timtot = 5000.d0
          timtot = 1100.d0
c         timtot = 1.50d0


c Setup tables for calculating exponentials
          call dexptablesmall_setup (dexptablesmall)
          call dexptablebig_setup   (dexptablebig)
          call otis_table_setup (otis_table,how_often,dt)

c Compartments contacted by "axoaxonic interneurons" are IS only
          do i = 1, num_L2pyr 
          do j = 1, num_placeholder2_to_L2pyr 
             com_placeholder2_to_L2pyr (j,i) = 69
          end do
          end do
          do i = 1, num_LOT
          do j = 1, num_placeholder2_to_LOT
             com_placeholder2_to_LOT(j,i) = 54
          end do
          end do
          do i = 1, num_semilunar   
          do j = 1, num_placeholder2_to_semilunar   
             com_placeholder2_to_semilunar   (j,i) = 69
          end do
          end do
c         do i = 1, num_multipolar   
c         do j = 1, num_placeholder2_to_multipolar   
c            com_placeholder2_to_multipolar   (j,i) = 56
c         end do
c         end do
          do i = 1, num_L3pyr   
          do j = 1, num_placeholder2_to_L3pyr   
             com_placeholder2_to_L3pyr   (j,i) = 69
          end do
          end do
! NOTE deepLTS was "deepaxax" in earlier code, so deepLTS
! connections need to be defined above
c         do i = 1, num_L2pyr    
c         do j = 1, num_deepLTS_to_L2pyr    
c            com_deepLTS_to_L2pyr    (j,i) = 69
c         end do
c         end do
c         do i = 1, num_LOT   
c         do j = 1, num_deepLTS_to_LOT   
c            com_deepLTS_to_LOT   (j,i) = 54
c         end do
c         end do
c         do i = 1, num_semilunar      
c         do j = 1, num_deepLTS_to_semilunar      
c            com_deepLTS_to_semilunar      (j,i) = 56 
c         end do
c         end do
c         do i = 1, num_multipolar      
c         do j = 1, num_deepLTS_to_multipolar      
c            com_deepLTS_to_multipolar      (j,i) = 56 
c         end do
c         end do
c         do i = 1, num_L3pyr      
c         do j = 1, num_deepLTS_to_L3pyr      
c            com_deepLTS_to_L3pyr      (j,i) = 45 
c         end do
c         end do
c End section on making axoaxonic cells connect to IS's

c Construct synaptic connectivity tables
                display = 0

          CALL synaptic_map_construct (thisno,
     &     num_L2pyr, num_L2pyr,           
     &     map_L2pyr_to_L2pyr,
     &     num_L2pyr_to_L2pyr,    display)
          CALL synaptic_map_construct (thisno,
     &     num_L2pyr, num_placeholder1,            
     &     map_L2pyr_to_placeholder1,  
     &     num_L2pyr_to_placeholder1,     display)
          CALL synaptic_map_construct (thisno,
     &     num_L2pyr, num_supng  ,            
     &     map_L2pyr_to_supng  ,  
     &     num_L2pyr_to_supng  ,     display)
          CALL synaptic_map_construct (thisno,
     &     num_L2pyr, num_placeholder2,            
     &     map_L2pyr_to_placeholder2,  
     &     num_L2pyr_to_placeholder2,     display)
          CALL synaptic_map_construct (thisno,
     &     num_L2pyr, num_placeholder3,             
     &     map_L2pyr_to_placeholder3,   
     &     num_L2pyr_to_placeholder3,      display)
          CALL synaptic_map_construct (thisno,
     &     num_L2pyr, num_LOT,          
     &     map_L2pyr_to_LOT,
     &     num_L2pyr_to_LOT,   display)
          CALL synaptic_map_construct (thisno,
     &     num_L2pyr, num_semilunar,             
     &     map_L2pyr_to_semilunar,   
     &     num_L2pyr_to_semilunar,      display)
          CALL synaptic_map_construct (thisno,
     &     num_L2pyr, num_multipolar,             
     &     map_L2pyr_to_multipolar,   
     &     num_L2pyr_to_multipolar,      display)
          CALL synaptic_map_construct (thisno,
     &     num_L2pyr, num_deepbask,           
     &     map_L2pyr_to_deepbask, 
     &     num_L2pyr_to_deepbask,    display)
          CALL synaptic_map_construct (thisno,
     &     num_L2pyr, num_deepLTS,           
     &     map_L2pyr_to_deepLTS, 
     &     num_L2pyr_to_deepLTS,    display)
          CALL synaptic_map_construct (thisno,
     &     num_L2pyr ,  num_deepng      ,             
     &     map_L2pyr_to_deepng      ,   
     &     num_L2pyr_to_deepng      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_L2pyr, num_supVIP ,            
     &     map_L2pyr_to_supVIP ,  
     &     num_L2pyr_to_supVIP ,     display)
          CALL synaptic_map_construct (thisno,
     &     num_L2pyr, num_L3pyr,          
     &     map_L2pyr_to_L3pyr,
     &     num_L2pyr_to_L3pyr,   display)

          CALL synaptic_map_construct (thisno,
     &     num_placeholder1, num_L2pyr,           
     &     map_placeholder1_to_L2pyr, 
     &     num_placeholder1_to_L2pyr,   display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder1, num_placeholder1,            
     &     map_placeholder1_to_placeholder1,  
     &     num_placeholder1_to_placeholder1,    display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder1, num_supng  ,            
     &     map_placeholder1_to_supng  ,  
     &     num_placeholder1_to_supng  ,    display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder1, num_placeholder2,            
     &     map_placeholder1_to_placeholder2,  
     &     num_placeholder1_to_placeholder2,    display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder1, num_placeholder3,             
     &     map_placeholder1_to_placeholder3,   
     &     num_placeholder1_to_placeholder3,     display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder1, num_LOT,          
     &     map_placeholder1_to_LOT,
     &     num_placeholder1_to_LOT,  display)

          CALL synaptic_map_construct (thisno,
     &     num_supng  , num_L2pyr ,          
     &     map_supng_to_L2pyr ,
     &     num_supng_to_L2pyr ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supng  , num_L3pyr,          
     &     map_supng_to_L3pyr,
     &     num_supng_to_L3pyr,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supng  , num_semilunar   ,          
     &     map_supng_to_semilunar   ,
     &     num_supng_to_semilunar   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supng  , num_multipolar   ,          
     &     map_supng_to_multipolar   ,
     &     num_supng_to_multipolar   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supng  , num_supng    ,          
     &     map_supng_to_supng    ,
     &     num_supng_to_supng    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supng  , num_placeholder1  ,          
     &     map_supng_to_placeholder1  ,
     &     num_supng_to_placeholder1  ,  display)

          CALL synaptic_map_construct (thisno,
     &     num_placeholder2, num_L2pyr,           
     &     map_placeholder2_to_L2pyr, 
     &     num_placeholder2_to_L2pyr,   display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder2, num_LOT,          
     &     map_placeholder2_to_LOT,
     &     num_placeholder2_to_LOT,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder2, num_semilunar,             
     &     map_placeholder2_to_semilunar,   
     &     num_placeholder2_to_semilunar,     display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder2, num_multipolar,             
     &     map_placeholder2_to_multipolar,   
     &     num_placeholder2_to_multipolar,     display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder2, num_L3pyr,             
     &     map_placeholder2_to_L3pyr,   
     &     num_placeholder2_to_L3pyr,  display)

          CALL synaptic_map_construct (thisno,
     &     num_placeholder3,  num_L2pyr,              
     &     map_placeholder3_to_L2pyr,    
     &     num_placeholder3_to_L2pyr ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder3,  num_placeholder1,               
     &     map_placeholder3_to_placeholder1,     
     &     num_placeholder3_to_placeholder1,    display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder3,  num_placeholder2,               
     &     map_placeholder3_to_placeholder2,     
     &     num_placeholder3_to_placeholder2,    display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder3,  num_placeholder3,                
     &     map_placeholder3_to_placeholder3,      
     &     num_placeholder3_to_placeholder3,     display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder3,  num_LOT,             
     &     map_placeholder3_to_LOT,   
     &     num_placeholder3_to_LOT,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder3,  num_semilunar,                
     &     map_placeholder3_to_semilunar,      
     &     num_placeholder3_to_semilunar   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder3,  num_multipolar,                
     &     map_placeholder3_to_multipolar,      
     &     num_placeholder3_to_multipolar   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder3,  num_deepbask,              
     &     map_placeholder3_to_deepbask,    
     &     num_placeholder3_to_deepbask ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder3,  num_deepLTS,              
     &     map_placeholder3_to_deepLTS,    
     &     num_placeholder3_to_deepLTS ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder3,  num_supVIP ,               
     &     map_placeholder3_to_supVIP ,     
     &     num_placeholder3_to_supVIP ,    display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder3,  num_L3pyr,             
     &     map_placeholder3_to_L3pyr,   
     &     num_placeholder3_to_L3pyr,  display)

          CALL synaptic_map_construct (thisno,
     &     num_LOT,  num_L2pyr,              
     &     map_LOT_to_L2pyr,    
     &     num_LOT_to_L2pyr,   display)
          CALL synaptic_map_construct (thisno,
     &     num_LOT,  num_placeholder1,               
     &     map_LOT_to_placeholder1,     
     &     num_LOT_to_placeholder1,    display)
          CALL synaptic_map_construct (thisno,
     &     num_LOT,  num_placeholder2,               
     &     map_LOT_to_placeholder2,     
     &     num_LOT_to_placeholder2,    display)
          CALL synaptic_map_construct (thisno,
     &     num_LOT,  num_placeholder3,                
     &     map_LOT_to_placeholder3,      
     &     num_LOT_to_placeholder3,     display)
          CALL synaptic_map_construct (thisno,
     &     num_LOT,  num_LOT,             
     &     map_LOT_to_LOT,   
     &     num_LOT_to_LOT,  display)
          CALL synaptic_map_construct (thisno,
     &     num_LOT,  num_semilunar,                
     &     map_LOT_to_semilunar,      
     &     num_LOT_to_semilunar,     display)
          CALL synaptic_map_construct (thisno,
     &     num_LOT,  num_multipolar,                
     &     map_LOT_to_multipolar,      
     &     num_LOT_to_multipolar,     display)
          CALL synaptic_map_construct (thisno,
     &     num_LOT,  num_deepbask,              
     &     map_LOT_to_deepbask,    
     &     num_LOT_to_deepbask,   display)
          CALL synaptic_map_construct (thisno,
     &     num_LOT,  num_deepng  ,              
     &     map_LOT_to_deepng  ,    
     &     num_LOT_to_deepng  ,   display)
          CALL synaptic_map_construct (thisno,
     &     num_LOT,  num_deepLTS,              
     &     map_LOT_to_deepLTS,    
     &     num_LOT_to_deepLTS,   display)
          CALL synaptic_map_construct (thisno,
     &     num_LOT,  num_supVIP ,               
     &     map_LOT_to_supVIP ,     
     &     num_LOT_to_supVIP ,    display)
          CALL synaptic_map_construct (thisno,
     &     num_LOT,  num_supng  ,               
     &     map_LOT_to_supng  ,     
     &     num_LOT_to_supng  ,    display)
          CALL synaptic_map_construct (thisno,
     &     num_LOT,  num_L3pyr,             
     &     map_LOT_to_L3pyr,   
     &     num_LOT_to_L3pyr,  display)

          CALL synaptic_map_construct (thisno,
     &     num_semilunar,  num_L2pyr,              
     &     map_semilunar_to_L2pyr,    
     &     num_semilunar_to_L2pyr,   display)
          CALL synaptic_map_construct (thisno,
     &     num_semilunar,  num_placeholder1,               
     &     map_semilunar_to_placeholder1,     
     &     num_semilunar_to_placeholder1,    display)
          CALL synaptic_map_construct (thisno,
     &     num_semilunar,  num_placeholder2,               
     &     map_semilunar_to_placeholder2,     
     &     num_semilunar_to_placeholder2,    display)
          CALL synaptic_map_construct (thisno,
     &     num_semilunar,  num_placeholder3,                
     &     map_semilunar_to_placeholder3,      
     &     num_semilunar_to_placeholder3,     display)
          CALL synaptic_map_construct (thisno,
     &     num_semilunar,  num_LOT,             
     &     map_semilunar_to_LOT,   
     &     num_semilunar_to_LOT,  display)
          CALL synaptic_map_construct (thisno,
     &     num_semilunar,  num_semilunar   ,             
     &     map_semilunar_to_semilunar   ,   
     &     num_semilunar_to_semilunar   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_semilunar,  num_multipolar   ,             
     &     map_semilunar_to_multipolar   ,   
     &     num_semilunar_to_multipolar   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_semilunar,  num_deepbask ,             
     &     map_semilunar_to_deepbask ,   
     &     num_semilunar_to_deepbask ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_semilunar,  num_deepng   ,             
     &     map_semilunar_to_deepng   ,   
     &     num_semilunar_to_deepng   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_semilunar,  num_deepLTS ,             
     &     map_semilunar_to_deepLTS ,   
     &     num_semilunar_to_deepLTS ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_semilunar,  num_supVIP   ,             
     &     map_semilunar_to_supVIP   ,   
     &     num_semilunar_to_supVIP   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_semilunar,  num_L3pyr,             
     &     map_semilunar_to_L3pyr,   
     &     num_semilunar_to_L3pyr,  display)

          CALL synaptic_map_construct (thisno,
     &     num_multipolar,  num_L2pyr ,             
     &     map_multipolar_to_L2pyr ,   
     &     num_multipolar_to_L2pyr ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_multipolar,  num_placeholder1  ,             
     &     map_multipolar_to_placeholder1  ,   
     &     num_multipolar_to_placeholder1  ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_multipolar,  num_placeholder2  ,             
     &     map_multipolar_to_placeholder2  ,   
     &     num_multipolar_to_placeholder2  ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_multipolar,  num_placeholder3   ,             
     &     map_multipolar_to_placeholder3   ,   
     &     num_multipolar_to_placeholder3   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_multipolar,  num_LOT,             
     &     map_multipolar_to_LOT,   
     &     num_multipolar_to_LOT,  display)
          CALL synaptic_map_construct (thisno,
     &     num_multipolar,  num_semilunar   ,             
     &     map_multipolar_to_semilunar   ,   
     &     num_multipolar_to_semilunar   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_multipolar,  num_multipolar   ,             
     &     map_multipolar_to_multipolar   ,   
     &     num_multipolar_to_multipolar   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_multipolar,  num_deepbask ,             
     &     map_multipolar_to_deepbask ,   
     &     num_multipolar_to_deepbask ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_multipolar,  num_deepng   ,             
     &     map_multipolar_to_deepng   ,   
     &     num_multipolar_to_deepng   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_multipolar,  num_deepLTS ,             
     &     map_multipolar_to_deepLTS ,   
     &     num_multipolar_to_deepLTS ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_multipolar,  num_supVIP   ,             
     &     map_multipolar_to_supVIP   ,   
     &     num_multipolar_to_supVIP   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_multipolar,  num_L3pyr,             
     &     map_multipolar_to_L3pyr,   
     &     num_multipolar_to_L3pyr,  display)

          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_LOT,             
     &     map_deepbask_to_LOT,   
     &     num_deepbask_to_LOT,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_semilunar   ,             
     &     map_deepbask_to_semilunar   ,   
     &     num_deepbask_to_semilunar   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_multipolar   ,             
     &     map_deepbask_to_multipolar   ,   
     &     num_deepbask_to_multipolar   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_deepbask ,             
     &     map_deepbask_to_deepbask ,   
     &     num_deepbask_to_deepbask ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_deepng   ,             
     &     map_deepbask_to_deepng   ,   
     &     num_deepbask_to_deepng   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_deepLTS ,             
     &     map_deepbask_to_deepLTS ,   
     &     num_deepbask_to_deepLTS ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_supVIP   ,             
     &     map_deepbask_to_supVIP   ,   
     &     num_deepbask_to_supVIP   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_L2pyr,             
     &     map_deepbask_to_L2pyr,   
     &     num_deepbask_to_L2pyr,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepbask,  num_L3pyr,             
     &     map_deepbask_to_L3pyr,   
     &     num_deepbask_to_L3pyr,  display)

          CALL synaptic_map_construct (thisno,
     &     num_deepng  ,  num_semilunar   ,             
     &     map_deepng_to_semilunar   ,   
     &     num_deepng_to_semilunar   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepng  ,  num_multipolar   ,             
     &     map_deepng_to_multipolar   ,   
     &     num_deepng_to_multipolar   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepng  ,  num_L2pyr,             
     &     map_deepng_to_L2pyr,   
     &     num_deepng_to_L2pyr,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepng  ,  num_L3pyr,             
     &     map_deepng_to_L3pyr,   
     &     num_deepng_to_L3pyr,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepng  ,  num_LOT,             
     &     map_deepng_to_LOT,   
     &     num_deepng_to_LOT,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepng  ,  num_deepng   ,             
     &     map_deepng_to_deepng   ,   
     &     num_deepng_to_deepng   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepng  ,  num_deepbask ,             
     &     map_deepng_to_deepbask ,   
     &     num_deepng_to_deepbask ,  display)

          CALL synaptic_map_construct (thisno,
     &     num_deepLTS,  num_L2pyr ,             
     &     map_deepLTS_to_L2pyr ,   
     &     num_deepLTS_to_L2pyr ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS,  num_LOT,             
     &     map_deepLTS_to_LOT,   
     &     num_deepLTS_to_LOT,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS,  num_semilunar   ,             
     &     map_deepLTS_to_semilunar   ,   
     &     num_deepLTS_to_semilunar   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS,  num_multipolar   ,             
     &     map_deepLTS_to_multipolar   ,   
     &     num_deepLTS_to_multipolar   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_deepLTS,  num_L3pyr   ,             
     &     map_deepLTS_to_L3pyr   ,   
     &     num_deepLTS_to_L3pyr   ,  display)

          CALL synaptic_map_construct (thisno,
     &     num_supVIP ,  num_L2pyr    ,             
     &     map_supVIP_to_L2pyr    ,   
     &     num_supVIP_to_L2pyr    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supVIP ,  num_placeholder1     ,             
     &     map_supVIP_to_placeholder1     ,   
     &     num_supVIP_to_placeholder1     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supVIP ,  num_placeholder2     ,             
     &     map_supVIP_to_placeholder2     ,   
     &     num_supVIP_to_placeholder2     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supVIP ,  num_placeholder3      ,             
     &     map_supVIP_to_placeholder3      ,   
     &     num_supVIP_to_placeholder3      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supVIP ,  num_supng       ,             
     &     map_supVIP_to_supng       ,   
     &     num_supVIP_to_supng       ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supVIP ,  num_LOT   ,             
     &     map_supVIP_to_LOT   ,   
     &     num_supVIP_to_LOT   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supVIP ,  num_semilunar      ,             
     &     map_supVIP_to_semilunar      ,   
     &     num_supVIP_to_semilunar      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supVIP ,  num_multipolar      ,             
     &     map_supVIP_to_multipolar      ,   
     &     num_supVIP_to_multipolar      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supVIP ,  num_deepbask    ,             
     &     map_supVIP_to_deepbask    ,   
     &     num_supVIP_to_deepbask    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supVIP ,  num_deepLTS    ,             
     &     map_supVIP_to_deepLTS    ,   
     &     num_supVIP_to_deepLTS    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supVIP ,  num_supVIP     ,             
     &     map_supVIP_to_supVIP     ,   
     &     num_supVIP_to_supVIP     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_supVIP ,  num_L3pyr   ,             
     &     map_supVIP_to_L3pyr   ,   
     &     num_supVIP_to_L3pyr   ,  display)

          CALL synaptic_map_construct (thisno,
     &     num_placeholder5 ,  num_L2pyr    ,             
     &     map_placeholder5_to_L2pyr    ,   
     &     num_placeholder5_to_L2pyr    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder5 ,  num_placeholder1     ,             
     &     map_placeholder5_to_placeholder1     ,   
     &     num_placeholder5_to_placeholder1     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder5 ,  num_supng       ,             
     &     map_placeholder5_to_supng       ,   
     &     num_placeholder5_to_supng       ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder5 ,  num_placeholder2     ,             
     &     map_placeholder5_to_placeholder2     ,   
     &     num_placeholder5_to_placeholder2     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder5 ,  num_LOT   ,             
     &     map_placeholder5_to_LOT   ,   
     &     num_placeholder5_to_LOT   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder5 ,  num_semilunar      ,             
     &     map_placeholder5_to_semilunar      ,   
     &     num_placeholder5_to_semilunar      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder5 ,  num_multipolar      ,             
     &     map_placeholder5_to_multipolar      ,   
     &     num_placeholder5_to_multipolar      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder5 ,  num_deepbask    ,             
     &     map_placeholder5_to_deepbask    ,   
     &     num_placeholder5_to_deepbask    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder5 ,  num_deepng      ,             
     &     map_placeholder5_to_deepng      ,   
     &     num_placeholder5_to_deepng      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder5 ,  num_deepLTS    ,             
     &     map_placeholder5_to_deepLTS    ,   
     &     num_placeholder5_to_deepLTS    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder5 ,  num_placeholder6         ,             
     &     map_placeholder5_to_placeholder6         ,   
     &     num_placeholder5_to_placeholder6         ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder5 ,  num_L3pyr   ,             
     &     map_placeholder5_to_L3pyr   ,   
     &     num_placeholder5_to_L3pyr   ,  display)

          CALL synaptic_map_construct (thisno,
     &     num_placeholder6 ,  num_placeholder5         ,             
     &     map_placeholder6_to_placeholder5         ,   
     &     num_placeholder6_to_placeholder5         ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_placeholder6 ,  num_placeholder6         ,             
     &     map_placeholder6_to_placeholder6         ,   
     &     num_placeholder6_to_placeholder6         ,  display)

          CALL synaptic_map_construct (thisno,
     &     num_L3pyr ,  num_L2pyr    ,             
     &     map_L3pyr_to_L2pyr    ,   
     &     num_L3pyr_to_L2pyr    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_L3pyr ,  num_placeholder1     ,             
     &     map_L3pyr_to_placeholder1     ,   
     &     num_L3pyr_to_placeholder1     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_L3pyr ,  num_placeholder2     ,             
     &     map_L3pyr_to_placeholder2     ,   
     &     num_L3pyr_to_placeholder2     ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_L3pyr ,  num_placeholder3      ,             
     &     map_L3pyr_to_placeholder3      ,   
     &     num_L3pyr_to_placeholder3      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_L3pyr ,  num_LOT   ,             
     &     map_L3pyr_to_LOT   ,   
     &     num_L3pyr_to_LOT   ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_L3pyr ,  num_semilunar      ,             
     &     map_L3pyr_to_semilunar      ,   
     &     num_L3pyr_to_semilunar      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_L3pyr ,  num_multipolar      ,             
     &     map_L3pyr_to_multipolar      ,   
     &     num_L3pyr_to_multipolar      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_L3pyr ,  num_deepbask    ,             
     &     map_L3pyr_to_deepbask    ,   
     &     num_L3pyr_to_deepbask    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_L3pyr ,  num_deepng      ,             
     &     map_L3pyr_to_deepng      ,   
     &     num_L3pyr_to_deepng      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_L3pyr ,  num_deepLTS    ,             
     &     map_L3pyr_to_deepLTS    ,   
     &     num_L3pyr_to_deepLTS    ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_L3pyr ,  num_supVIP      ,             
     &     map_L3pyr_to_supVIP      ,   
     &     num_L3pyr_to_supVIP      ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_L3pyr ,  num_placeholder5         ,             
     &     map_L3pyr_to_placeholder5         ,   
     &     num_L3pyr_to_placeholder5         ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_L3pyr ,  num_placeholder6         ,             
     &     map_L3pyr_to_placeholder6         ,   
     &     num_L3pyr_to_placeholder6         ,  display)
          CALL synaptic_map_construct (thisno,
     &     num_L3pyr ,  num_L3pyr   ,             
     &     map_L3pyr_to_L3pyr   ,   
     &     num_L3pyr_to_L3pyr   ,  display)
c Finish construction of synaptic connection tables.

c Count and display how many to-be activated LOT inputs each cell has
             if (thisno.eq.0) then
           write (6,1497)
1497       format('  number activated inputs to L2pyr cells')
           do L = 1, num_L2pyr
            ictr = 0
            do k = 1, num_LOT_to_L2pyr
             j = map_LOT_to_L2pyr (k,L)
             if (j.le.25) ictr = ictr + 1
            end do
             if (ictr.gt.0) then
            write (6,1498) L, ictr
             endif
1498        format(2i7)
           end do

           write (6,1499)
1499       format('  number activated inputs to L3pyr cells')
           do L = 1, num_L3pyr
            ictr = 0
            do k = 1, num_LOT_to_L3pyr
             j = map_LOT_to_L3pyr (k,L)
             if (j.le.25) ictr = ictr + 1
            end do
             if (ictr.gt.0) then
            write (6,1498) L, ictr
             endif
           end do

           write (6,1492)
1492       format('  number activated inputs to SL cells')
           do L = 1, num_semilunar
            ictr = 0
            do k = 1, num_LOT_to_semilunar
             j = map_LOT_to_semilunar (k,L)
             if (j.le.25) ictr = ictr + 1
            end do
             if (ictr.gt.0) then
            write (6,1498) L, ictr
             endif
           end do
             endif ! thisno = 0

c Construct synaptic compartment maps.  
                display = 0

          CALL synaptic_compmap_construct (thisno,
     &     num_L2pyr, com_L2pyr_to_L2pyr,           
     &     num_L2pyr_to_L2pyr,
     &  ncompallow_L2pyr_to_L2pyr,
     &   compallow_L2pyr_to_L2pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder1  , com_L2pyr_to_placeholder1,            
     &     num_L2pyr_to_placeholder1,
     &    ncompallow_L2pyr_to_placeholder1,  
     &     compallow_L2pyr_to_placeholder1,   display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supng    , com_L2pyr_to_supng  ,            
     &     num_L2pyr_to_supng  ,
     &    ncompallow_L2pyr_to_supng  ,  
     &     compallow_L2pyr_to_supng  ,   display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder2  , com_L2pyr_to_placeholder2,            
     &     num_L2pyr_to_placeholder2,  
     &    ncompallow_L2pyr_to_placeholder2,  
     &     compallow_L2pyr_to_placeholder2,   display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder3   , com_L2pyr_to_placeholder3,             
     &     num_L2pyr_to_placeholder3,   
     &    ncompallow_L2pyr_to_placeholder3,   
     &     compallow_L2pyr_to_placeholder3,    display)

          CALL synaptic_compmap_construct (thisno,
     &     num_LOT, com_L2pyr_to_LOT,          
     &     num_L2pyr_to_LOT,
     &    ncompallow_L2pyr_to_LOT,
     &     compallow_L2pyr_to_LOT, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_semilunar   , com_L2pyr_to_semilunar   ,          
     &     num_L2pyr_to_semilunar   ,
     &    ncompallow_L2pyr_to_semilunar   ,
     &     compallow_L2pyr_to_semilunar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_multipolar   , com_L2pyr_to_multipolar   ,          
     &     num_L2pyr_to_multipolar   ,
     &    ncompallow_L2pyr_to_multipolar   ,
     &     compallow_L2pyr_to_multipolar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_L2pyr_to_deepbask ,          
     &     num_L2pyr_to_deepbask ,
     &    ncompallow_L2pyr_to_deepbask ,
     &     compallow_L2pyr_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS , com_L2pyr_to_deepLTS ,          
     &     num_L2pyr_to_deepLTS ,
     &    ncompallow_L2pyr_to_deepLTS ,
     &     compallow_L2pyr_to_deepLTS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_L2pyr_to_deepng   ,          
     &     num_L2pyr_to_deepng   ,
     &    ncompallow_L2pyr_to_deepng   ,
     &     compallow_L2pyr_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supVIP  , com_L2pyr_to_supVIP  ,          
     &     num_L2pyr_to_supVIP  ,
     &    ncompallow_L2pyr_to_supVIP  ,
     &     compallow_L2pyr_to_supVIP  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L3pyr, com_L2pyr_to_L3pyr,          
     &     num_L2pyr_to_L3pyr,
     &    ncompallow_L2pyr_to_L3pyr,
     &     compallow_L2pyr_to_L3pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L2pyr , com_placeholder1_to_L2pyr ,          
     &     num_placeholder1_to_L2pyr ,
     &    ncompallow_placeholder1_to_L2pyr ,
     &     compallow_placeholder1_to_L2pyr , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder1  , com_placeholder1_to_placeholder1  ,          
     &     num_placeholder1_to_placeholder1  ,
     &    ncompallow_placeholder1_to_placeholder1  ,
     &     compallow_placeholder1_to_placeholder1  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supng    , com_placeholder1_to_supng    ,          
     &     num_placeholder1_to_supng    ,
     &    ncompallow_placeholder1_to_supng    ,
     &     compallow_placeholder1_to_supng    , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder2  , com_placeholder1_to_placeholder2  ,          
     &     num_placeholder1_to_placeholder2  ,
     &    ncompallow_placeholder1_to_placeholder2  ,
     &     compallow_placeholder1_to_placeholder2  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder3   , com_placeholder1_to_placeholder3   ,          
     &     num_placeholder1_to_placeholder3   ,
     &    ncompallow_placeholder1_to_placeholder3   ,
     &     compallow_placeholder1_to_placeholder3   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_LOT, com_placeholder1_to_LOT,          
     &     num_placeholder1_to_LOT,
     &    ncompallow_placeholder1_to_LOT,
     &     compallow_placeholder1_to_LOT, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L2pyr , com_supng_to_L2pyr ,          
     &     num_supng_to_L2pyr ,
     &    ncompallow_supng_to_L2pyr ,
     &     compallow_supng_to_L2pyr , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L3pyr, com_supng_to_L3pyr,          
     &     num_supng_to_L3pyr,
     &    ncompallow_supng_to_L3pyr,
     &     compallow_supng_to_L3pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_semilunar   , com_supng_to_semilunar   ,          
     &     num_supng_to_semilunar   ,
     &    ncompallow_supng_to_semilunar   ,
     &     compallow_supng_to_semilunar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_multipolar   , com_supng_to_multipolar   ,          
     &     num_supng_to_multipolar   ,
     &    ncompallow_supng_to_multipolar   ,
     &     compallow_supng_to_multipolar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supng    , com_supng_to_supng    ,          
     &     num_supng_to_supng    ,
     &    ncompallow_supng_to_supng    ,
     &     compallow_supng_to_supng    , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder1  , com_supng_to_placeholder1  ,          
     &     num_supng_to_placeholder1  ,
     &    ncompallow_supng_to_placeholder1  ,
     &     compallow_supng_to_placeholder1  , display)

! Calls below not necessary?
          CALL synaptic_compmap_construct (thisno,
     &     num_L2pyr , com_placeholder2_to_L2pyr ,          
     &     num_placeholder2_to_L2pyr ,
     &    ncompallow_placeholder2_to_L2pyr ,
     &     compallow_placeholder2_to_L2pyr , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_LOT, com_placeholder2_to_LOT,          
     &     num_placeholder2_to_LOT,
     &    ncompallow_placeholder2_to_LOT,
     &     compallow_placeholder2_to_LOT, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_semilunar   , com_placeholder2_to_semilunar   ,          
     &     num_placeholder2_to_semilunar   ,
     &    ncompallow_placeholder2_to_semilunar   ,
     &     compallow_placeholder2_to_semilunar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_multipolar   , com_placeholder2_to_multipolar   ,          
     &     num_placeholder2_to_multipolar   ,
     &    ncompallow_placeholder2_to_multipolar   ,
     &     compallow_placeholder2_to_multipolar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L3pyr, com_placeholder2_to_L3pyr,          
     &     num_placeholder2_to_L3pyr,
     &    ncompallow_placeholder2_to_L3pyr,
     &     compallow_placeholder2_to_L3pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L2pyr , com_placeholder3_to_L2pyr ,          
     &     num_placeholder3_to_L2pyr ,
     &    ncompallow_placeholder3_to_L2pyr ,
     &     compallow_placeholder3_to_L2pyr , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder1  , com_placeholder3_to_placeholder1  ,          
     &     num_placeholder3_to_placeholder1  ,
     &    ncompallow_placeholder3_to_placeholder1  ,
     &     compallow_placeholder3_to_placeholder1  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder2  , com_placeholder3_to_placeholder2  ,          
     &     num_placeholder3_to_placeholder2  ,
     &    ncompallow_placeholder3_to_placeholder2  ,
     &     compallow_placeholder3_to_placeholder2  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder3   , com_placeholder3_to_placeholder3   ,          
     &     num_placeholder3_to_placeholder3   ,
     &    ncompallow_placeholder3_to_placeholder3   ,
     &     compallow_placeholder3_to_placeholder3   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_LOT, com_placeholder3_to_LOT,          
     &     num_placeholder3_to_LOT,
     &    ncompallow_placeholder3_to_LOT,
     &     compallow_placeholder3_to_LOT, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_semilunar   , com_placeholder3_to_semilunar   ,          
     &     num_placeholder3_to_semilunar   ,
     &    ncompallow_placeholder3_to_semilunar   ,
     &     compallow_placeholder3_to_semilunar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_multipolar   , com_placeholder3_to_multipolar   ,          
     &     num_placeholder3_to_multipolar   ,
     &    ncompallow_placeholder3_to_multipolar   ,
     &     compallow_placeholder3_to_multipolar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_placeholder3_to_deepbask ,          
     &     num_placeholder3_to_deepbask ,
     &    ncompallow_placeholder3_to_deepbask ,
     &     compallow_placeholder3_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS , com_placeholder3_to_deepLTS ,          
     &     num_placeholder3_to_deepLTS ,
     &    ncompallow_placeholder3_to_deepLTS ,
     &     compallow_placeholder3_to_deepLTS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supVIP  , com_placeholder3_to_supVIP  ,          
     &     num_placeholder3_to_supVIP  ,
     &    ncompallow_placeholder3_to_supVIP  ,
     &     compallow_placeholder3_to_supVIP  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L3pyr, com_placeholder3_to_L3pyr,          
     &     num_placeholder3_to_L3pyr,
     &    ncompallow_placeholder3_to_L3pyr,
     &     compallow_placeholder3_to_L3pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L2pyr , com_LOT_to_L2pyr ,          
     &     num_LOT_to_L2pyr ,
     &    ncompallow_LOT_to_L2pyr ,
     &     compallow_LOT_to_L2pyr , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder1  , com_LOT_to_placeholder1  ,          
     &     num_LOT_to_placeholder1  ,
     &    ncompallow_LOT_to_placeholder1  ,
     &     compallow_LOT_to_placeholder1  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder2  , com_LOT_to_placeholder2  ,          
     &     num_LOT_to_placeholder2  ,
     &    ncompallow_LOT_to_placeholder2  ,
     &     compallow_LOT_to_placeholder2  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder3   , com_LOT_to_placeholder3   ,          
     &     num_LOT_to_placeholder3   ,
     &    ncompallow_LOT_to_placeholder3   ,
     &     compallow_LOT_to_placeholder3   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_LOT, com_LOT_to_LOT,          
     &     num_LOT_to_LOT,
     &    ncompallow_LOT_to_LOT,
     &     compallow_LOT_to_LOT, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_semilunar   , com_LOT_to_semilunar   ,          
     &     num_LOT_to_semilunar   ,
     &    ncompallow_LOT_to_semilunar   ,
     &     compallow_LOT_to_semilunar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_multipolar   , com_LOT_to_multipolar   ,          
     &     num_LOT_to_multipolar   ,
     &    ncompallow_LOT_to_multipolar   ,
     &     compallow_LOT_to_multipolar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_LOT_to_deepbask ,          
     &     num_LOT_to_deepbask ,
     &    ncompallow_LOT_to_deepbask ,
     &     compallow_LOT_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_LOT_to_deepng   ,          
     &     num_LOT_to_deepng   ,
     &    ncompallow_LOT_to_deepng   ,
     &     compallow_LOT_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS , com_LOT_to_deepLTS ,          
     &     num_LOT_to_deepLTS ,
     &    ncompallow_LOT_to_deepLTS ,
     &     compallow_LOT_to_deepLTS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supVIP  , com_LOT_to_supVIP  ,          
     &     num_LOT_to_supVIP  ,
     &    ncompallow_LOT_to_supVIP  ,
     &     compallow_LOT_to_supVIP  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supng   , com_LOT_to_supng   ,          
     &     num_LOT_to_supng   ,
     &    ncompallow_LOT_to_supng   ,
     &     compallow_LOT_to_supng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L3pyr, com_LOT_to_L3pyr,          
     &     num_LOT_to_L3pyr,
     &    ncompallow_LOT_to_L3pyr,
     &     compallow_LOT_to_L3pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L2pyr , com_semilunar_to_L2pyr ,          
     &     num_semilunar_to_L2pyr ,
     &    ncompallow_semilunar_to_L2pyr ,
     &     compallow_semilunar_to_L2pyr , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder1  , com_semilunar_to_placeholder1  ,          
     &     num_semilunar_to_placeholder1  ,
     &    ncompallow_semilunar_to_placeholder1  ,
     &     compallow_semilunar_to_placeholder1  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder2  , com_semilunar_to_placeholder2  ,          
     &     num_semilunar_to_placeholder2  ,
     &    ncompallow_semilunar_to_placeholder2  ,
     &     compallow_semilunar_to_placeholder2  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder3   , com_semilunar_to_placeholder3   ,          
     &     num_semilunar_to_placeholder3   ,
     &    ncompallow_semilunar_to_placeholder3   ,
     &     compallow_semilunar_to_placeholder3   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_LOT, com_semilunar_to_LOT,          
     &     num_semilunar_to_LOT,
     &    ncompallow_semilunar_to_LOT,
     &     compallow_semilunar_to_LOT, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_semilunar   , com_semilunar_to_semilunar   ,          
     &     num_semilunar_to_semilunar   ,
     &    ncompallow_semilunar_to_semilunar   ,
     &     compallow_semilunar_to_semilunar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_multipolar   , com_semilunar_to_multipolar   ,          
     &     num_semilunar_to_multipolar   ,
     &    ncompallow_semilunar_to_multipolar   ,
     &     compallow_semilunar_to_multipolar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_semilunar_to_deepbask ,          
     &     num_semilunar_to_deepbask ,
     &    ncompallow_semilunar_to_deepbask ,
     &     compallow_semilunar_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_semilunar_to_deepng   ,          
     &     num_semilunar_to_deepng   ,
     &    ncompallow_semilunar_to_deepng   ,
     &     compallow_semilunar_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS , com_semilunar_to_deepLTS ,          
     &     num_semilunar_to_deepLTS ,
     &    ncompallow_semilunar_to_deepLTS ,
     &     compallow_semilunar_to_deepLTS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supVIP  , com_semilunar_to_supVIP  ,          
     &     num_semilunar_to_supVIP  ,
     &    ncompallow_semilunar_to_supVIP  ,
     &     compallow_semilunar_to_supVIP  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L3pyr, com_semilunar_to_L3pyr,          
     &     num_semilunar_to_L3pyr,
     &    ncompallow_semilunar_to_L3pyr,
     &     compallow_semilunar_to_L3pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L2pyr , com_multipolar_to_L2pyr ,          
     &     num_multipolar_to_L2pyr ,
     &    ncompallow_multipolar_to_L2pyr ,
     &     compallow_multipolar_to_L2pyr , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder1  , com_multipolar_to_placeholder1  ,          
     &     num_multipolar_to_placeholder1  ,
     &    ncompallow_multipolar_to_placeholder1  ,
     &     compallow_multipolar_to_placeholder1  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder2  , com_multipolar_to_placeholder2  ,          
     &     num_multipolar_to_placeholder2  ,
     &    ncompallow_multipolar_to_placeholder2  ,
     &     compallow_multipolar_to_placeholder2  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder3   , com_multipolar_to_placeholder3   ,          
     &     num_multipolar_to_placeholder3   ,
     &    ncompallow_multipolar_to_placeholder3   ,
     &     compallow_multipolar_to_placeholder3   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_LOT, com_multipolar_to_LOT,          
     &     num_multipolar_to_LOT,
     &    ncompallow_multipolar_to_LOT,
     &     compallow_multipolar_to_LOT, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_semilunar   , com_multipolar_to_semilunar   ,          
     &     num_multipolar_to_semilunar   ,
     &    ncompallow_multipolar_to_semilunar   ,
     &     compallow_multipolar_to_semilunar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_multipolar   , com_multipolar_to_multipolar   ,          
     &     num_multipolar_to_multipolar   ,
     &    ncompallow_multipolar_to_multipolar   ,
     &     compallow_multipolar_to_multipolar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_multipolar_to_deepbask ,          
     &     num_multipolar_to_deepbask ,
     &    ncompallow_multipolar_to_deepbask ,
     &     compallow_multipolar_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_multipolar_to_deepng   ,          
     &     num_multipolar_to_deepng   ,
     &    ncompallow_multipolar_to_deepng   ,
     &     compallow_multipolar_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS , com_multipolar_to_deepLTS ,          
     &     num_multipolar_to_deepLTS ,
     &    ncompallow_multipolar_to_deepLTS ,
     &     compallow_multipolar_to_deepLTS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supVIP  , com_multipolar_to_supVIP  ,          
     &     num_multipolar_to_supVIP  ,
     &    ncompallow_multipolar_to_supVIP  ,
     &     compallow_multipolar_to_supVIP  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L3pyr, com_multipolar_to_L3pyr,          
     &     num_multipolar_to_L3pyr,
     &    ncompallow_multipolar_to_L3pyr,
     &     compallow_multipolar_to_L3pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_LOT, com_deepbask_to_LOT,          
     &     num_deepbask_to_LOT,
     &    ncompallow_deepbask_to_LOT,
     &     compallow_deepbask_to_LOT, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_semilunar   , com_deepbask_to_semilunar   ,          
     &     num_deepbask_to_semilunar   ,
     &    ncompallow_deepbask_to_semilunar   ,
     &     compallow_deepbask_to_semilunar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_multipolar   , com_deepbask_to_multipolar   ,        
     &     num_deepbask_to_multipolar   ,
     &    ncompallow_deepbask_to_multipolar   ,
     &     compallow_deepbask_to_multipolar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_deepbask_to_deepbask ,          
     &     num_deepbask_to_deepbask ,
     &    ncompallow_deepbask_to_deepbask ,
     &     compallow_deepbask_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_deepbask_to_deepng   ,          
     &     num_deepbask_to_deepng   ,
     &    ncompallow_deepbask_to_deepng   ,
     &     compallow_deepbask_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS , com_deepbask_to_deepLTS ,          
     &     num_deepbask_to_deepLTS ,
     &    ncompallow_deepbask_to_deepLTS ,
     &     compallow_deepbask_to_deepLTS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supVIP  , com_deepbask_to_supVIP  ,          
     &     num_deepbask_to_supVIP  ,
     &    ncompallow_deepbask_to_supVIP  ,
     &     compallow_deepbask_to_supVIP  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L2pyr, com_deepbask_to_L2pyr,          
     &     num_deepbask_to_L2pyr,
     &    ncompallow_deepbask_to_L2pyr,
     &     compallow_deepbask_to_L2pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L3pyr, com_deepbask_to_L3pyr,          
     &     num_deepbask_to_L3pyr,
     &    ncompallow_deepbask_to_L3pyr,
     &     compallow_deepbask_to_L3pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_semilunar   , com_deepng_to_semilunar   ,          
     &     num_deepng_to_semilunar   ,
     &    ncompallow_deepng_to_semilunar   ,
     &     compallow_deepng_to_semilunar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_multipolar   , com_deepng_to_multipolar   ,          
     &     num_deepng_to_multipolar   ,
     &    ncompallow_deepng_to_multipolar   ,
     &     compallow_deepng_to_multipolar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L2pyr, com_deepng_to_L2pyr,          
     &     num_deepng_to_L2pyr,
     &    ncompallow_deepng_to_L2pyr,
     &     compallow_deepng_to_L2pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L3pyr, com_deepng_to_L3pyr,          
     &     num_deepng_to_L3pyr,
     &    ncompallow_deepng_to_L3pyr,
     &     compallow_deepng_to_L3pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_LOT, com_deepng_to_LOT,          
     &     num_deepng_to_LOT,
     &    ncompallow_deepng_to_LOT,
     &     compallow_deepng_to_LOT, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_deepng_to_deepng   ,          
     &     num_deepng_to_deepng   ,
     &    ncompallow_deepng_to_deepng   ,
     &     compallow_deepng_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_deepng_to_deepbask ,          
     &     num_deepng_to_deepbask ,
     &    ncompallow_deepng_to_deepbask ,
     &     compallow_deepng_to_deepbask , display)

! Below calls not necessary in plateau, but now ARE necessary
          CALL synaptic_compmap_construct (thisno,
     &     num_L2pyr , com_deepLTS_to_L2pyr ,          
     &     num_deepLTS_to_L2pyr ,
     &    ncompallow_deepLTS_to_L2pyr ,
     &     compallow_deepLTS_to_L2pyr , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_LOT, com_deepLTS_to_LOT,          
     &     num_deepLTS_to_LOT,
     &    ncompallow_deepLTS_to_LOT,
     &     compallow_deepLTS_to_LOT, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_semilunar   , com_deepLTS_to_semilunar   ,          
     &     num_deepLTS_to_semilunar   ,
     &    ncompallow_deepLTS_to_semilunar   ,
     &     compallow_deepLTS_to_semilunar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_multipolar   , com_deepLTS_to_multipolar   ,          
     &     num_deepLTS_to_multipolar   ,
     &    ncompallow_deepLTS_to_multipolar   ,
     &     compallow_deepLTS_to_multipolar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L3pyr, com_deepLTS_to_L3pyr,          
     &     num_deepLTS_to_L3pyr,
     &    ncompallow_deepLTS_to_L3pyr,
     &     compallow_deepLTS_to_L3pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L2pyr , com_supVIP_to_L2pyr ,          
     &     num_supVIP_to_L2pyr ,
     &    ncompallow_supVIP_to_L2pyr ,
     &     compallow_supVIP_to_L2pyr , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder1  , com_supVIP_to_placeholder1  ,          
     &     num_supVIP_to_placeholder1  ,
     &    ncompallow_supVIP_to_placeholder1  ,
     &     compallow_supVIP_to_placeholder1  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder2  , com_supVIP_to_placeholder2  ,          
     &     num_supVIP_to_placeholder2  ,
     &    ncompallow_supVIP_to_placeholder2  ,
     &     compallow_supVIP_to_placeholder2  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder3   , com_supVIP_to_placeholder3   ,          
     &     num_supVIP_to_placeholder3   ,
     &    ncompallow_supVIP_to_placeholder3   ,
     &     compallow_supVIP_to_placeholder3   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supng    , com_supVIP_to_supng    ,          
     &     num_supVIP_to_supng    ,
     &    ncompallow_supVIP_to_supng    ,
     &     compallow_supVIP_to_supng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_LOT, com_supVIP_to_LOT,          
     &     num_supVIP_to_LOT,
     &    ncompallow_supVIP_to_LOT,
     &     compallow_supVIP_to_LOT, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_semilunar   , com_supVIP_to_semilunar   ,          
     &     num_supVIP_to_semilunar   ,
     &    ncompallow_supVIP_to_semilunar   ,
     &     compallow_supVIP_to_semilunar   , display)
!  Make connections explicit, so one cell to one compartment
c         do L = 1, num_semilunar
c         do i = 1, num_supVIP_to_semilunar ! should equal ncompallow
c          com_supVIP_to_semilunar (i,L) = 
c    &        compallow_supVIP_to_semilunar (i)
c         end do
c         end do

          CALL synaptic_compmap_construct (thisno,
     &     num_multipolar   , com_supVIP_to_multipolar   ,          
     &     num_supVIP_to_multipolar   ,
     &    ncompallow_supVIP_to_multipolar   ,
     &     compallow_supVIP_to_multipolar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_supVIP_to_deepbask ,          
     &     num_supVIP_to_deepbask ,
     &    ncompallow_supVIP_to_deepbask ,
     &     compallow_supVIP_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS , com_supVIP_to_deepLTS ,          
     &     num_supVIP_to_deepLTS ,
     &    ncompallow_supVIP_to_deepLTS ,
     &     compallow_supVIP_to_deepLTS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supVIP  , com_supVIP_to_supVIP  ,          
     &     num_supVIP_to_supVIP  ,
     &    ncompallow_supVIP_to_supVIP  ,
     &     compallow_supVIP_to_supVIP  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L3pyr, com_supVIP_to_L3pyr,          
     &     num_supVIP_to_L3pyr,
     &    ncompallow_supVIP_to_L3pyr,
     &     compallow_supVIP_to_L3pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L2pyr , com_placeholder5_to_L2pyr ,          
     &     num_placeholder5_to_L2pyr ,
     &    ncompallow_placeholder5_to_L2pyr ,
     &     compallow_placeholder5_to_L2pyr , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder1  , com_placeholder5_to_placeholder1  ,          
     &     num_placeholder5_to_placeholder1  ,
     &    ncompallow_placeholder5_to_placeholder1  ,
     &     compallow_placeholder5_to_placeholder1  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supng    , com_placeholder5_to_supng    ,          
     &     num_placeholder5_to_supng    ,
     &    ncompallow_placeholder5_to_supng    ,
     &     compallow_placeholder5_to_supng    , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder2  , com_placeholder5_to_placeholder2  ,          
     &     num_placeholder5_to_placeholder2  ,
     &    ncompallow_placeholder5_to_placeholder2  ,
     &     compallow_placeholder5_to_placeholder2  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_LOT, com_placeholder5_to_LOT,          
     &     num_placeholder5_to_LOT,
     &    ncompallow_placeholder5_to_LOT,
     &     compallow_placeholder5_to_LOT, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_semilunar   , com_placeholder5_to_semilunar   ,          
     &     num_placeholder5_to_semilunar   ,
     &    ncompallow_placeholder5_to_semilunar   ,
     &     compallow_placeholder5_to_semilunar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_multipolar   , com_placeholder5_to_multipolar   ,          
     &     num_placeholder5_to_multipolar   ,
     &    ncompallow_placeholder5_to_multipolar   ,
     &     compallow_placeholder5_to_multipolar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_placeholder5_to_deepbask ,          
     &     num_placeholder5_to_deepbask ,
     &    ncompallow_placeholder5_to_deepbask ,
     &     compallow_placeholder5_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_placeholder5_to_deepng   ,          
     &     num_placeholder5_to_deepng   ,
     &    ncompallow_placeholder5_to_deepng   ,
     &     compallow_placeholder5_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS , com_placeholder5_to_deepLTS ,          
     &     num_placeholder5_to_deepLTS ,
     &    ncompallow_placeholder5_to_deepLTS ,
     &     compallow_placeholder5_to_deepLTS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder6      , com_placeholder5_to_placeholder6   ,          
     &     num_placeholder5_to_placeholder6      ,
     &    ncompallow_placeholder5_to_placeholder6      ,
     &     compallow_placeholder5_to_placeholder6      , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L3pyr, com_placeholder5_to_L3pyr,          
     &     num_placeholder5_to_L3pyr,
     &    ncompallow_placeholder5_to_L3pyr,
     &     compallow_placeholder5_to_L3pyr, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder5      , com_placeholder6_to_placeholder5   ,          
     &     num_placeholder6_to_placeholder5      ,
     &    ncompallow_placeholder6_to_placeholder5      ,
     &     compallow_placeholder6_to_placeholder5      , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder6      , com_placeholder6_to_placeholder6   ,          
     &     num_placeholder6_to_placeholder6      ,
     &    ncompallow_placeholder6_to_placeholder6      ,
     &     compallow_placeholder6_to_placeholder6      , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L2pyr , com_L3pyr_to_L2pyr ,          
     &     num_L3pyr_to_L2pyr ,
     &    ncompallow_L3pyr_to_L2pyr ,
     &     compallow_L3pyr_to_L2pyr , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder1  , com_L3pyr_to_placeholder1  ,          
     &     num_L3pyr_to_placeholder1  ,
     &    ncompallow_L3pyr_to_placeholder1  ,
     &     compallow_L3pyr_to_placeholder1  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder2  , com_L3pyr_to_placeholder2  ,          
     &     num_L3pyr_to_placeholder2  ,
     &    ncompallow_L3pyr_to_placeholder2  ,
     &     compallow_L3pyr_to_placeholder2  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder3   , com_L3pyr_to_placeholder3   ,          
     &     num_L3pyr_to_placeholder3   ,
     &    ncompallow_L3pyr_to_placeholder3   ,
     &     compallow_L3pyr_to_placeholder3   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_LOT, com_L3pyr_to_LOT,          
     &     num_L3pyr_to_LOT,
     &    ncompallow_L3pyr_to_LOT,
     &     compallow_L3pyr_to_LOT, display)

          CALL synaptic_compmap_construct (thisno,
     &     num_semilunar   , com_L3pyr_to_semilunar   ,          
     &     num_L3pyr_to_semilunar   ,
     &    ncompallow_L3pyr_to_semilunar   ,
     &     compallow_L3pyr_to_semilunar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_multipolar   , com_L3pyr_to_multipolar   ,          
     &     num_L3pyr_to_multipolar   ,
     &    ncompallow_L3pyr_to_multipolar   ,
     &     compallow_L3pyr_to_multipolar   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepbask , com_L3pyr_to_deepbask ,          
     &     num_L3pyr_to_deepbask ,
     &    ncompallow_L3pyr_to_deepbask ,
     &     compallow_L3pyr_to_deepbask , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepng   , com_L3pyr_to_deepng   ,          
     &     num_L3pyr_to_deepng   ,
     &    ncompallow_L3pyr_to_deepng   ,
     &     compallow_L3pyr_to_deepng   , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_deepLTS , com_L3pyr_to_deepLTS ,          
     &     num_L3pyr_to_deepLTS ,
     &    ncompallow_L3pyr_to_deepLTS ,
     &     compallow_L3pyr_to_deepLTS , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_supVIP  , com_L3pyr_to_supVIP  ,          
     &     num_L3pyr_to_supVIP  ,
     &    ncompallow_L3pyr_to_supVIP  ,
     &     compallow_L3pyr_to_supVIP  , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder5      , com_L3pyr_to_placeholder5 ,          
     &     num_L3pyr_to_placeholder5      ,
     &    ncompallow_L3pyr_to_placeholder5      ,
     &     compallow_L3pyr_to_placeholder5      , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_placeholder6      , com_L3pyr_to_placeholder6      ,          
     &     num_L3pyr_to_placeholder6      ,
     &    ncompallow_L3pyr_to_placeholder6      ,
     &     compallow_L3pyr_to_placeholder6      , display)

          CALL synaptic_compmap_construct (thisno,
     &     num_L3pyr, com_L3pyr_to_L3pyr,          
     &     num_L3pyr_to_L3pyr,
     &    ncompallow_L3pyr_to_L3pyr,
     &     compallow_L3pyr_to_L3pyr, display)

c Finish construction of synaptic compartment maps. 


c Construct gap-junction tables
! axax interneurons a special case
c          gjtable_placeholder2(1,1) = 1
c          gjtable_placeholder2(1,2) = 12
c          gjtable_placeholder2(1,3) = 2
c          gjtable_placeholder2(1,4) = 12
           gjtable_deepLTS(1,1) = 1
           gjtable_deepLTS(1,2) = 12
           gjtable_deepLTS(1,3) = 2
           gjtable_deepLTS(1,4) = 12

c CALL BELOW WILL HAVE TO BE ADJUSTED, SO CONNECTIONS REMAIN WITHIN NODE 
c ASSUME 2 NODES FOR L2pyr
      CALL groucho_gapbld (thisno, num_L2pyr,
     & totaxgj_L2pyr  , gjtable_L2pyr,
     & table_axgjcompallow_L2pyr, 
     & num_axgjcompallow_L2pyr, 0) 
         do i = 1, totaxgj_L2pyr
           j = gjtable_L2pyr (i,1)
           k = gjtable_L2pyr (i,3)
       if ((j.le.ncellspernode).and.(k.gt.ncellspernode)) then
        gjtable_L2pyr(i,3) = gjtable_L2pyr(i,3) - ncellspernode 
       else if ((j.gt.ncellspernode).and.(k.le.ncellspernode)) then
        gjtable_L2pyr(i,3) = gjtable_L2pyr(i,3) + ncellspernode 
       endif
         end do  ! i

      CALL groucho_gapbld (thisno, num_LOT,
     & totaxgj_LOT , gjtable_LOT,
     & table_axgjcompallow_LOT,
     & num_axgjcompallow_LOT,  0) 

      CALL groucho_gapbld (thisno, num_semilunar,   
     & totaxgj_semilunar    , gjtable_semilunar   ,
     & table_axgjcompallow_semilunar   ,
     & num_axgjcompallow_semilunar   ,  0) 

      CALL groucho_gapbld (thisno, num_multipolar,   
     & totaxgj_multipolar    , gjtable_multipolar   ,
     & table_axgjcompallow_multipolar   ,
     & num_axgjcompallow_multipolar   ,  0) 

      CALL groucho_gapbld (thisno, num_L3pyr,   
     & totaxgj_L3pyr    , gjtable_L3pyr   ,
     & table_axgjcompallow_L3pyr   ,
     & num_axgjcompallow_L3pyr   ,  0) 

      CALL groucho_gapbld (thisno, num_placeholder1  ,   
     & totSDgj_placeholder1      , gjtable_placeholder1     ,
     & table_SDgjcompallow_placeholder1     ,
     & num_SDgjcompallow_placeholder1     ,  0) 

      CALL groucho_gapbld (thisno, num_supng    ,   
     & totSDgj_supng        , gjtable_supng       ,
     & table_SDgjcompallow_supng       ,
     & num_SDgjcompallow_supng       ,  0) 

      CALL groucho_gapbld (thisno, num_placeholder3   ,   
     & totSDgj_placeholder3       , gjtable_placeholder3      ,
     & table_SDgjcompallow_placeholder3      ,
     & num_SDgjcompallow_placeholder3      ,  0) 

      CALL groucho_gapbld (thisno, num_deepbask ,   
     & totSDgj_deepbask     , gjtable_deepbask    ,
     & table_SDgjcompallow_deepbask    ,
     & num_SDgjcompallow_deepbask    ,  0) 

      CALL groucho_gapbld (thisno, num_deepng   ,   
     & totSDgj_deepng       , gjtable_deepng      ,
     & table_SDgjcompallow_deepng      ,
     & num_SDgjcompallow_deepng      ,  0) 

      CALL groucho_gapbld (thisno, num_supVIP  ,   
     & totSDgj_supVIP      , gjtable_supVIP     ,
     & table_SDgjcompallow_supVIP     ,
     & num_SDgjcompallow_supVIP     ,  0) 

      CALL groucho_gapbld (thisno, num_placeholder5      ,   
     & totaxgj_placeholder5          , gjtable_placeholder5         ,
     & table_axgjcompallow_placeholder5         ,
     & num_axgjcompallow_placeholder5         ,  0) 

      CALL groucho_gapbld (thisno, num_placeholder6      ,   
     & totaxgj_placeholder6          , gjtable_placeholder6         ,
     & table_axgjcompallow_placeholder6         ,
     & num_axgjcompallow_placeholder6         ,  0) 


! Define tonic currents to different cell types
       call durand(seed,num_L2pyr ,ranvec_L2pyr )
       do L = 1, num_L2pyr 
       do i = 69, 74  ! axonal compartments
c       curr_L2pyr  (i,L) = -0.013d0 
        curr_L2pyr  (i,L) = -0.016d0 
c       curr_L2pyr  (i,L) =  0.050d0 + 0.005d0 *
          if (L.le.500) then
c       curr_L2pyr  (1,L) =  0.00d0 + 0.15d0 *
        curr_L2pyr  (1,L) =  0.20d0 + 0.15d0 *
     &      dble (L) / 500.d0  ! note gradient 
          else
c       curr_L2pyr  (1,L) =  0.00d0 + 0.15d0 *
        curr_L2pyr  (1,L) =  0.20d0 + 0.15d0 *
     &      dble (L-500) / 500.d0  ! note gradient 
          endif
c    &    ranvec_L2pyr (L)
       end do
       end do
c      curr_L2pyr (1,4) = 0.15d0  ! DEPOLARIZE 1 CELL

c      call durand(seed,num_placeholder1  ,ranvec_placeholder1  )
c      do L = 1, num_placeholder1    
c       curr_placeholder1   (1,L) = -0.10d0 + 0.02d0 *
c       curr_placeholder1   (1,L) = -0.04d0 + 0.02d0 *
c    &    ranvec_placeholder1  (L)
c      end do

c      call durand(seed,num_placeholder2  ,ranvec_placeholder2  )
c      do L = 1, num_placeholder2    
c       curr_placeholder2   (1,L) = -0.10d0 + 0.02d0 *
c       curr_placeholder2   (1,L) = -0.04d0 + 0.02d0 *
c    &    ranvec_placeholder2  (L)
c      end do

       call durand(seed,num_supVIP  ,ranvec_supVIP   )
       do L = 1, num_supVIP    
c       curr_supVIP    (1,L) = -0.10d0 + 0.02d0 *
c       curr_supVIP    (1,L) =  0.040d0 + 0.01d0 *
        curr_supVIP    (1,L) =  0.130d0 + 0.01d0 *
     &    ranvec_supVIP (L)
       end do

       call durand(seed,num_deepbask  ,ranvec_deepbask  )
       do L = 1, num_deepbask    
        curr_deepbask   (1,L) = -0.10d0 + 0.02d0 *
     &    ranvec_deepbask  (L)
       end do

       do L = 1, num_supng
          curr_supng (1,L) = -0.03d0 ! to suppress spontaneous firing
       end do

       call durand(seed,num_LOT,ranvec_LOT)
       do L = 1, num_LOT  
        curr_LOT (1,L) = -0.25d0 + 0.05d0 *
c       curr_LOT (1,L) =  0.00d0 + 0.05d0 *
     &    ranvec_LOT(L)
       end do

       call durand(seed,num_semilunar   ,ranvec_semilunar   )
       do L = 1, num_semilunar    
          do i = 2, 37
c       curr_semilunar    (i,L) = 0.055d0 + 0.001d0 *  ! current to basal/oblique dendrites 
c       curr_semilunar    (i,L) = 0.015d0 + 0.001d0 *  ! current to basal/oblique dendrites 
        curr_semilunar    (i,L) = 0.005d0 + 0.001d0 *  ! current to basal/oblique dendrites 
! but note that basal dendrites disconnected
     &    ranvec_semilunar   (L)
          end do
         do i = 41, 41  ! parts of shaft, tuft 
           curr_semilunar (i,L) = 0.05d0
         end do
c       curr_semilunar (57,L) = -0.015d0 ! axon
        curr_semilunar (70,L) =  0.000d0 ! axon
        curr_semilunar (71,L) =  0.000d0 ! axon
        curr_semilunar (72,L) =  0.000d0 ! axon
        curr_semilunar (73,L) =  0.000d0 ! axon
        curr_semilunar (74,L) =  0.000d0 ! axon
       end do

c      call durand(seed,num_multipolar   ,ranvec_multipolar   )
c      do L = 1, num_multipolar    
c       curr_multipolar    (1,L) = -0.1d0 + 0.02d0 * 
c    &    dble (L) / dble (num_multipolar)   ! note gradient        
c    &         ranvec_multipolar (L)
c         do i = 2, 34
c       curr_multipolar    (i,L) = -0.01d0 + 0.01d0 *  ! current to basal/oblique dendrites 
c       curr_multipolar    (i,L) =  0.00d0 + 0.01d0 *  ! current to basal/oblique dendrites 
c    &    ranvec_multipolar   (L)
c         end do
c       curr_multipolar (54,L) = -0.01d0 ! axon
c       curr_multipolar (55,L) = -0.01d0 ! axon
c       curr_multipolar (56,L) = -0.01d0 ! axon
c       curr_multipolar (57,L) = -0.01d0 ! axon
c       curr_multipolar (58,L) = -0.01d0 ! axon
c       curr_multipolar (59,L) = -0.01d0 ! axon
c      end do

       do L = 1, num_multipolar
        curr_multipolar (1,L) = -0.250d0
       end do

       do L = 1, num_supng
          curr_supng (1,L) = -0.03d0 ! to suppress spontaneous firing
       end do

       do L = 1, num_deepng
c         curr_deepng (1,L) = -0.025d0 ! to suppress spontaneous firing
          curr_deepng (1,L) = -0.045d0 ! to suppress spontaneous firing
       end do

       call durand(seed,num_L3pyr ,ranvec_L3pyr )
c            z = 0.70d0
c            z = 0.10d0
c            z = -0.10d0
             z = -0.20d0
       do L = 1, num_L3pyr
c       curr_L3pyr  (1,L) = z - 0.2d0 * ! decrease spont. firing
c       curr_L3pyr  (1,L) = z - 0.1d0 * ! decrease spont. firing
        curr_L3pyr  (1,L) = z + 0.1d0 * ! decrease spont. firing
     &    ranvec_L3pyr (L)
       end do

c      call durand(seed,num_placeholder6       ,ranvec_placeholder6   )
c      do L = 1, num_placeholder6          
c       curr_placeholder6        (1,L) = -0.05d0 + 0.05d0 *
c    &    ranvec_placeholder6       (L)
c      end do

c      call durand(seed,num_placeholder5       ,ranvec_placeholder5   )
c      do L = 1, num_placeholder5          
c       curr_placeholder5        (1,L) = 0.00d0 + 0.01d0 *
c    &    ranvec_placeholder5       (L)
c      end do

! ? remove from the picture some cells, by hyperpolarizing the respective axons
           go to 9901
             do L = 1, num_L2pyr
              curr_L2pyr(numcomp_L2pyr,L) = -0.25d0
             end do

             do L = 1, num_placeholder1  
              curr_placeholder1  (numcomp_placeholder1  ,L) = -0.25d0
             end do

             do L = 1, num_supng    
              curr_supng    (numcomp_supng    ,L) = -0.25d0
             end do

             do L = 1, num_placeholder2  
              curr_placeholder2  (numcomp_placeholder2  ,L) = -0.25d0
             end do

             do L = 1, num_placeholder3   
              curr_placeholder3   (numcomp_placeholder3   ,L) = -0.25d0
             end do

             do L = 1, num_supVIP   
              curr_supVIP   (numcomp_supVIP   ,L) = -0.25d0
             end do

             do L = 1, num_LOT
              curr_LOT(numcomp_LOT,L) = -0.25d0
             end do
9901       continue

! ? remove from the picture other cells, by hyperpolarizing the respective axons
           go to 9902
             do L = 1, num_semilunar   
              curr_semilunar   (numcomp_semilunar   ,L) = -0.25d0
             end do

             do L = 1, num_multipolar   
              curr_multipolar   (numcomp_multipolar   ,L) = -0.25d0
             end do

             do L = 1, num_L3pyr
              curr_L3pyr(numcomp_L3pyr,L) = -0.25d0
             end do

             do L = 1, num_deepbask 
              curr_deepbask (numcomp_deepbask ,L) = -0.25d0
             end do

             do L = 1, num_deepng   
              curr_deepng   (numcomp_deepng   ,L) = -0.25d0
             end do

             do L = 1, num_deepLTS 
              curr_deepLTS (numcomp_deepLTS ,L) = -0.25d0
             end do

             do L = 1, num_supVIP   
              curr_supVIP   (numcomp_supVIP   ,L) = -0.25d0
             end do
9902       continue

       seed = 137.d0

       O = 0
       time = 0.d0


c INITIALIZE ALL THE INTEGRATION SUBROUTINES
        initialize = 0
        firstcell = 1
        lastcell =  1
      IF (nodecell(thisno).eq.'L2pyr       ') then
       CALL INTEGRATE_L2pyr       (O, time, num_L2pyr,
     &    V_L2pyr, curr_L2pyr,
     &    initialize, firstcell, lastcell,
     & gAMPA_L2pyr, gNMDA_L2pyr, gGABA_A_L2pyr,
     & gGABA_B_L2pyr, Mg, 
     & gapcon_L2pyr  ,totaxgj_L2pyr   ,gjtable_L2pyr, dt,
     &  chi_L2pyr,mnaf_L2pyr,mnap_L2pyr,
     &  hnaf_L2pyr,mkdr_L2pyr,mka_L2pyr,
     &  hka_L2pyr,mk2_L2pyr,hk2_L2pyr,
     &  mkm_L2pyr,mkc_L2pyr,mkahp_L2pyr,
     &  mcat_L2pyr,hcat_L2pyr,mcal_L2pyr,
     &  mar_L2pyr,field_sup   ,field_deep,rel_axonshift_L2pyr)

c     ELSE if (nodecell(thisno).eq.'placeholder1  ') then
      ELSE if (nodecell(thisno).eq.'supintern   ') then
       CALL INTEGRATE_placeholder1 (O, time, num_placeholder1 ,
     &    V_placeholder1 , curr_placeholder1 ,
     $    initialize, firstcell, lastcell,
     & gAMPA_placeholder1 , gNMDA_placeholder1 , gGABA_A_placeholder1 ,
     & Mg, 
     & gapcon_placeholder1,totSDgj_placeholder1,gjtable_placeholder1,dt,
     &  chi_placeholder1,mnaf_placeholder1,mnap_placeholder1,
     &  hnaf_placeholder1,mkdr_placeholder1,mka_placeholder1,
     &  hka_placeholder1,mk2_placeholder1,hk2_placeholder1,
     &  mkm_placeholder1,mkc_placeholder1,mkahp_placeholder1,
     &  mcat_placeholder1,hcat_placeholder1,mcal_placeholder1,
     &  mar_placeholder1)


c     ELSE if (nodecell(thisno).eq.'supng    ') then
       CALL INTEGRATE_supng    (O, time, num_supng   ,
     &    V_supng   , curr_supng   ,
     $    initialize, firstcell, lastcell,
     & gAMPA_supng   , gNMDA_supng   , gGABA_A_supng   ,
     & Mg, 
     & gapcon_supng     ,totSDgj_supng      ,gjtable_supng   , dt,
     &  chi_supng  ,mnaf_supng  ,mnap_supng  ,
     &  hnaf_supng  ,mkdr_supng  ,mka_supng  ,
     &  hka_supng  ,mk2_supng  ,hk2_supng  ,
     &  mkm_supng  ,mkc_supng  ,mkahp_supng  ,
     &  mcat_supng  ,hcat_supng  ,mcal_supng  ,
     &  mar_supng  )

c     ELSE if (nodecell(thisno).eq.'placeholder2  ') then
       CALL INTEGRATE_placeholder2 (O, time, num_placeholder2 ,
     &    V_placeholder2 , curr_placeholder2 ,
     &    initialize, firstcell, lastcell,
     & gAMPA_placeholder2, gNMDA_placeholder2 , gGABA_A_placeholder2 ,
     & Mg, 
     & gapcon_placeholder2,totSDgj_placeholder2,gjtable_placeholder2,dt,
     &  chi_placeholder2,mnaf_placeholder2,mnap_placeholder2,
     &  hnaf_placeholder2,mkdr_placeholder2,mka_placeholder2,
     &  hka_placeholder2,mk2_placeholder2,hk2_placeholder2,
     &  mkm_placeholder2,mkc_placeholder2,mkahp_placeholder2,
     &  mcat_placeholder2,hcat_placeholder2,mcal_placeholder2,
     &  mar_placeholder2)

c     ELSE if (nodecell(thisno).eq.'placeholder3   ') then
       CALL INTEGRATE_placeholder3  (O, time, num_placeholder3  ,
     &    V_placeholder3  , curr_placeholder3  ,
     &    initialize, firstcell, lastcell,
     & gAMPA_placeholder3  , gNMDA_placeholder3  ,gGABA_A_placeholder3,
     & Mg, 
     & gapcon_placeholder3,totSDgj_placeholder3,gjtable_placeholder3,dt,
     &  chi_placeholder3,mnaf_placeholder3,mnap_placeholder3,
     &  hnaf_placeholder3,mkdr_placeholder3,mka_placeholder3,
     &  hka_placeholder3,mk2_placeholder3,hk2_placeholder3,
     &  mkm_placeholder3,mkc_placeholder3,mkahp_placeholder3,
     &  mcat_placeholder3,hcat_placeholder3,mcal_placeholder3,
     &  mar_placeholder3)

       CALL INTEGRATE_supVIP   (O, time, num_supVIP  ,
     &    V_supVIP  , curr_supVIP  ,
     & initialize, firstcell, lastcell,
     & gAMPA_supVIP  , gNMDA_supVIP  , gGABA_A_supVIP  ,
     & Mg, 
     & gapcon_supVIP  ,totSDgj_supVIP  ,gjtable_supVIP  , dt,
     &  chi_supVIP,mnaf_supVIP,mnap_supVIP,
     &  hnaf_supVIP,mkdr_supVIP,mka_supVIP,
     &  hka_supVIP,mk2_supVIP,hk2_supVIP,
     &  mkm_supVIP,mkc_supVIP,mkahp_supVIP,
     &  mcat_supVIP,hcat_supVIP,mcal_supVIP,
     &  mar_supVIP)

      ELSE if (nodecell(thisno).eq.'LOT         ') then
       CALL INTEGRATE_LOT              (O, time, num_LOT,
     &    V_LOT, curr_LOT,
     &    initialize, firstcell, lastcell,
     & gAMPA_LOT, gNMDA_LOT, gGABA_A_LOT,
     & gGABA_B_LOT, Mg, 
     & gapcon_LOT,totaxgj_LOT,gjtable_LOT, dt,
     &  chi_LOT,mnaf_LOT,mnap_LOT,
     &  hnaf_LOT,mkdr_LOT,mka_LOT,
     &  hka_LOT,mk2_LOT,hk2_LOT,
     &  mkm_LOT,mkc_LOT,mkahp_LOT,
     &  mcat_LOT,hcat_LOT,mcal_LOT,
     &  mar_LOT)

       ELSE IF (nodecell(thisno).eq.'semilunar   ') then
       CALL INTEGRATE_semilunar       (O, time, num_semilunar,
     &    V_semilunar, curr_semilunar,
     &    initialize, firstcell, lastcell,
     & gAMPA_semilunar, gNMDA_semilunar, gGABA_A_semilunar,
     & gGABA_B_semilunar, Mg, 
     & gapcon_semilunar  ,totaxgj_semilunar   ,gjtable_semilunar, dt,
     &  chi_semilunar,mnaf_semilunar,mnap_semilunar,
     &  hnaf_semilunar,mkdr_semilunar,mka_semilunar,
     &  hka_semilunar,mk2_semilunar,hk2_semilunar,
     &  mkm_semilunar,mkc_semilunar,mkahp_semilunar,
     &  mcat_semilunar,hcat_semilunar,mcal_semilunar,
     &  mar_semilunar,field_sup   ,field_deep,rel_axonshift_semilunar)

      ELSE if (nodecell(thisno).eq.'multipolar') then
       CALL INTEGRATE_multipolar (O, time, num_multipolar,
     &    V_multipolar, curr_multipolar,
     & initialize, firstcell, lastcell,
     & gAMPA_multipolar, gNMDA_multipolar, gGABA_A_multipolar,
c    & gGABA_B_multipolar, Mg, 
     &                     Mg, 
     & gapcon_multipolar,totaxgj_multipolar,gjtable_multipolar,dt,
     &  chi_multipolar,mnaf_multipolar,mnap_multipolar,
     &  hnaf_multipolar,mkdr_multipolar,mka_multipolar,
     &  hka_multipolar,mk2_multipolar,hk2_multipolar,
     &  mkm_multipolar,mkc_multipolar,mkahp_multipolar,
     &  mcat_multipolar,hcat_multipolar,mcal_multipolar,
c    &  mar_multipolar,field_sup       ,field_deep       )
     &  mar_multipolar       )

      ELSE IF (nodecell(thisno).eq.'L3pyr       ') then
       CALL INTEGRATE_L3pyr       (O, time, num_L3pyr,
     &    V_L3pyr, curr_L3pyr,
     &    initialize, firstcell, lastcell,
     & gAMPA_L3pyr, gNMDA_L3pyr, gGABA_A_L3pyr,
     & gGABA_B_L3pyr, Mg, 
     & gapcon_L3pyr  ,totaxgj_L3pyr   ,gjtable_L3pyr, dt,
     &  chi_L3pyr,mnaf_L3pyr,mnap_L3pyr,
     &  hnaf_L3pyr,mkdr_L3pyr,mka_L3pyr,
     &  hka_L3pyr,mk2_L3pyr,hk2_L3pyr,
     &  mkm_L3pyr,mkc_L3pyr,mkahp_L3pyr,
     &  mcat_L3pyr,hcat_L3pyr,mcal_L3pyr,
     &  mar_L3pyr,field_sup   ,field_deep,rel_axonshift_L3pyr)

c     ELSE if (nodecell(thisno).eq.'deepbask ') then
      ELSE if (nodecell(thisno).eq.'deepintern  ') then
       CALL INTEGRATE_deepbask  (O, time, num_deepbask ,
     &    V_deepbask , curr_deepbask ,
     & initialize, firstcell, lastcell,
     & gAMPA_deepbask, gNMDA_deepbask, gGABA_A_deepbask,
     & Mg, 
     & gapcon_deepbask  ,totSDgj_deepbask   ,gjtable_deepbask, dt,
     &  chi_deepbask,mnaf_deepbask,mnap_deepbask,
     &  hnaf_deepbask,mkdr_deepbask,mka_deepbask,
     &  hka_deepbask,mk2_deepbask,hk2_deepbask,
     &  mkm_deepbask,mkc_deepbask,mkahp_deepbask,
     &  mcat_deepbask,hcat_deepbask,mcal_deepbask,
     &  mar_deepbask)

c     ELSE if (nodecell(thisno).eq.'deepng   ') then
       CALL INTEGRATE_deepng     (O, time, num_deepng   ,
     &    V_deepng   , curr_deepng   ,
     & initialize, firstcell, lastcell,
     & gAMPA_deepng  , gNMDA_deepng  , gGABA_A_deepng  ,
     & Mg, 
     & gapcon_deepng    ,totSDgj_deepng     ,gjtable_deepng  , dt,
     &  chi_deepng  ,mnaf_deepng  ,mnap_deepng  ,
     &  hnaf_deepng  ,mkdr_deepng  ,mka_deepng  ,
     &  hka_deepng  ,mk2_deepng  ,hk2_deepng  ,
     &  mkm_deepng  ,mkc_deepng  ,mkahp_deepng  ,
     &  mcat_deepng  ,hcat_deepng  ,mcal_deepng  ,
     &  mar_deepng  )

c     ELSE if (nodecell(thisno).eq.'deepLTS ') then
       CALL INTEGRATE_deepLTS   (O, time, num_deepLTS ,
     &    V_deepLTS , curr_deepLTS ,
     & initialize, firstcell, lastcell,
     & gAMPA_deepLTS, gNMDA_deepLTS, gGABA_A_deepLTS,
     & Mg, 
     & gapcon_deepLTS  ,totSDgj_deepLTS   ,gjtable_deepLTS, dt,
     &  chi_deepLTS,mnaf_deepLTS,mnap_deepLTS,
     &  hnaf_deepLTS,mkdr_deepLTS,mka_deepLTS,
     &  hka_deepLTS,mk2_deepLTS,hk2_deepLTS,
     &  mkm_deepLTS,mkc_deepLTS,mkahp_deepLTS,
     &  mcat_deepLTS,hcat_deepLTS,mcal_deepLTS,
     &  mar_deepLTS)

      ELSE if (nodecell(thisno).eq.'placeholder5') then
       CALL INTEGRATE_placeholder5      (O, time, num_placeholder5    ,
     &    V_placeholder5      , curr_placeholder5      ,
     & initialize, firstcell, lastcell,
     & gAMPA_placeholder5,gNMDA_placeholder5, gGABA_A_placeholder5    ,
     & gGABA_B_placeholder5, Mg, 
     & gapcon_placeholder5,totaxgj_placeholder5,gjtable_placeholder5,dt,
     &  chi_placeholder5,mnaf_placeholder5,mnap_placeholder5,
     &  hnaf_placeholder5,mkdr_placeholder5,mka_placeholder5,
     &  hka_placeholder5,mk2_placeholder5,hk2_placeholder5,
     &  mkm_placeholder5,mkc_placeholder5,mkahp_placeholder5,
     &  mcat_placeholder5,hcat_placeholder5,mcal_placeholder5,
     &  mar_placeholder5)

      ELSE if (nodecell(thisno).eq.'placeholder6') then
       CALL INTEGRATE_placeholder6 (O, time, num_placeholder6      ,
     &    V_placeholder6      , curr_placeholder6      ,
     & initialize, firstcell, lastcell,
     & gAMPA_placeholder6,gNMDA_placeholder6,gGABA_A_placeholder6   ,
     & gGABA_B_placeholder6, Mg, 
     & gapcon_placeholder6,totaxgj_placeholder6,gjtable_placeholder6,dt,
     &  chi_placeholder6,mnaf_placeholder6,mnap_placeholder6,
     &  hnaf_placeholder6,mkdr_placeholder6,mka_placeholder6,
     &  hka_placeholder6,mk2_placeholder6,hk2_placeholder6,
     &  mkm_placeholder6,mkc_placeholder6,mkahp_placeholder6,
     &  mcat_placeholder6,hcat_placeholder6,mcal_placeholder6,
     &  mar_placeholder6,field_sup,field_deep,rel_axonshift_L2pyr)

      ENDIF
c END INITIALIZATION OF INTEGRATION SUBROUTINES
c Note how superficial and deep interneuron calls lumped together


c BEGIN guts of main program.
c Each node takes care of all the cells of a particular type.
c On a node: enumerate the cells of its type; calculate their
c  synaptic inputs; set applied currents, including those
c  required by ectopic generation; call the numerical integration
c  subroutine; set up the distal_axon vector.  Each node 
c  broadcasts its own distal_axon vector to all the others, and also
c  receives distal_axon vectors from all the others.
c Then, update outtime array and outctr vector.  Repeat.

c CHANGE SEED?
        seed = 293.d0

1000    O = O + 1
        time = time + dt
        if (time.gt.timtot) goto 2000

c         gapcon_L2pyr = 12.00d-3
          gapcon_L2pyr =  0.d-3
c         gapcon_L2pyr =  5.d-3

c ? current pulse to some cells ?
c         do L = 1, num_L2pyr
c       if ((time.gt.600.d0).and.(time.lt.1100.d0)) then
c        curr_L2pyr(1,L) = 1.0d0 + dble(L) * 0.5d0/1000.d0
c       else
c        curr_L2pyr(1,L) = 0.d0
c       endif
c         end do

c current pulse to some multipolar
c       if ((time.gt.600.d0).and.(time.le.605.d0)) then
c         do L = 1, 20
c          curr_multipolar (1,L) = 0.50d0
c         end do
c         do L = 21, num_multipolar
c          curr_multipolar (1,L) = -0.1
c         end do
c       else
c         do L = 1, num_multipolar
c          curr_multipolar (1,L) = -0.1d0
c         end do
c       endif


c Define shift of semilunar axonal gNa rate functions, & other axon shifts
       rel_axonshift_semilunar = 5.0d0 + 0.0d0 * time/timtot
       rel_axonshift_L2pyr = 5.d0                      
       rel_axonshift_L3pyr = 5.d0                      

       noisepe_semilunar = noisepe_semilunar_save
       noisepe_multipolar = noisepe_multipolar_save

       initialize = 1  ! so integration subroutines actually integrate


c      IF (THISNO.EQ.0) THEN
       IF (nodecell(thisno) .eq. 'L2pyr       ') THEN
c L2pyr

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 + (i-1) * ncellspernode
          lastcell = firstcell - 1 + ncellspernode

          IF (MOD(O,how_often).eq.0) then
c 1st set L2pyr synaptic conductances to 0:

          do i = 1, numcomp_L2pyr
!         do j = 1, num_L2pyr
          do j = firstcell, lastcell ! Note
         gAMPA_L2pyr(i,j)   = 0.d0
         gNMDA_L2pyr(i,j)   = 0.d0
         gGABA_A_L2pyr(i,j) = 0.d0
         gGABA_B_L2pyr(i,j) = 0.d0
          end do
          end do

!        do L = 1, num_L2pyr
         do L = firstcell, lastcell  ! Note

c Handle L2pyr   -> L2pyr
      do i = 1, num_L2pyr_to_L2pyr
       j = map_L2pyr_to_L2pyr(i,L) ! j = presynaptic cell
       k = com_L2pyr_to_L2pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L2pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L2pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L2pyr_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_L2pyr(k,L)  = gAMPA_L2pyr(k,L) +
     &  gAMPA_L2pyr_to_L2pyr * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_L2pyr(k,L) = gNMDA_L2pyr(k,L) +
     &  gNMDA_L2pyr_to_L2pyr * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L2pyr_to_L2pyr
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_L2pyr(k,L) = gNMDA_L2pyr(k,L) +
     &  gNMDA_L2pyr_to_L2pyr * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L2pyr_to_L2pyr
       if (gNMDA_L2pyr(k,L).gt.z)
     &  gNMDA_L2pyr(k,L) = z
! end NMDA part

       end do ! m
      end do ! i



c Handle placeholder1    -> L2pyr
      do i = 1, num_placeholder1_to_L2pyr
       j = map_placeholder1_to_L2pyr(i,L) ! j = presynaptic cell
       k = com_placeholder1_to_L2pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder1(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder1(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder1_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L2pyr(k,L)  = gGABA_A_L2pyr(k,L) +
     &  gGABA_placeholder1_to_L2pyr * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle supng      -> L2pyr
      do i = 1, num_supng_to_L2pyr
       j = map_supng_to_L2pyr(i,L) ! j = presynaptic cell
       k = com_supng_to_L2pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supng(j)  ! enumerate presyn. spikes
        presyntime = outtime_supng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part ! NOTE DIFFERENT GABA-B HERE, NOW
        dexparg = delta / tauGABA_supng_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L2pyr(k,L)  = gGABA_A_L2pyr(k,L) +
     &  gGABA_supng_to_L2pyr * z      
! end GABA-A part

        dexparg = delta / tauGABAB_supng_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_B_L2pyr(k,L)  = gGABA_A_L2pyr(k,L) +
     &  gGABAB_supng_to_L2pyr * z      
! end GABA-A part


! end GABA-B part

       end do ! m
      end do ! i

c Handle deepbask   -> L2pyr
      do i = 1, num_deepbask_to_L2pyr   
       j = map_deepbask_to_L2pyr(i,L) ! j = presynaptic cell
       k = com_deepbask_to_L2pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_L2pyr   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L2pyr(k,L)  = gGABA_A_L2pyr(k,L) +
     &  gGABA_deepbask_to_L2pyr * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepng      -> L2pyr
      do i = 1, num_deepng_to_L2pyr
       j = map_deepng_to_L2pyr(i,L) ! j = presynaptic cell
       k = com_deepng_to_L2pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepng(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_deepng_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L2pyr(k,L)  = gGABA_A_L2pyr(k,L) +
     &  gGABA_deepng_to_L2pyr * z      
! end GABA-A part

        dexparg = delta / tauGABAB_deepng_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_B_L2pyr(k,L)  = gGABA_B_L2pyr(k,L) +
     &  gGABAB_deepng_to_L2pyr * z      

! end GABA-B part

       end do ! m
      end do ! i



c Handle placeholder2    -> L2pyr
      do i = 1, num_placeholder2_to_L2pyr
       j = map_placeholder2_to_L2pyr(i,L) ! j = presynaptic cell
       k = com_placeholder2_to_L2pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder2(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder2(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder2_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L2pyr(k,L)  = gGABA_A_L2pyr(k,L) +
     &  gGABA_placeholder2_to_L2pyr * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle placeholder3     -> L2pyr
      do i = 1, num_placeholder3_to_L2pyr
       j = map_placeholder3_to_L2pyr(i,L) ! j = presynaptic cell
       k = com_placeholder3_to_L2pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder3(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder3(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder3_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L2pyr(k,L)  = gGABA_A_L2pyr(k,L) +
     &  gGABA_placeholder3_to_L2pyr * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle LOT  -> L2pyr
      do i = 1, num_LOT_to_L2pyr
       j = map_LOT_to_L2pyr(i,L) ! j = presynaptic cell
       k = com_LOT_to_L2pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_LOT(j)  ! enumerate presyn. spikes
        presyntime = outtime_LOT(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_LOT_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_L2pyr(k,L)  = gAMPA_L2pyr(k,L) +
     &  gAMPA_LOT_to_L2pyr * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_L2pyr(k,L) = gNMDA_L2pyr(k,L) +
     &  gNMDA_LOT_to_L2pyr * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_LOT_to_L2pyr
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_L2pyr(k,L) = gNMDA_L2pyr(k,L) +
     &  gNMDA_LOT_to_L2pyr * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_LOT_to_L2pyr
       if (gNMDA_L2pyr(k,L).gt.z)
     &  gNMDA_L2pyr(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle semilunar     -> L2pyr
      do i = 1, num_semilunar_to_L2pyr
       j = map_semilunar_to_L2pyr(i,L) ! j = presynaptic cell
       k = com_semilunar_to_L2pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_semilunar(j)  ! enumerate presyn. spikes
        presyntime = outtime_semilunar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_semilunar_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_L2pyr(k,L)  = gAMPA_L2pyr(k,L) +
     &  gAMPA_semilunar_to_L2pyr * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_L2pyr(k,L) = gNMDA_L2pyr(k,L) +
     &  gNMDA_semilunar_to_L2pyr * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_semilunar_to_L2pyr
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_L2pyr(k,L) = gNMDA_L2pyr(k,L) +
     &  gNMDA_semilunar_to_L2pyr * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_semilunar_to_L2pyr
       if (gNMDA_L2pyr(k,L).gt.z)
     &  gNMDA_L2pyr(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle multipolar       -> L2pyr
      do i = 1, num_multipolar_to_L2pyr
       j = map_multipolar_to_L2pyr(i,L) ! j = presynaptic cell
       k = com_multipolar_to_L2pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_multipolar(j)  ! enumerate presyn. spikes
        presyntime = outtime_multipolar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_multipolar_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_L2pyr(k,L)  = gAMPA_L2pyr(k,L) +
     &  gAMPA_multipolar_to_L2pyr * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_L2pyr(k,L) = gNMDA_L2pyr(k,L) +
     &  gNMDA_multipolar_to_L2pyr * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_multipolar_to_L2pyr
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_L2pyr(k,L) = gNMDA_L2pyr(k,L) +
     &  gNMDA_multipolar_to_L2pyr * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_multipolar_to_L2pyr
       if (gNMDA_L2pyr(k,L).gt.z)
     &  gNMDA_L2pyr(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepLTS   -> L2pyr
      do i = 1, num_deepLTS_to_L2pyr
       j = map_deepLTS_to_L2pyr(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_L2pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L2pyr(k,L)  = gGABA_A_L2pyr(k,L) +
     &  gGABA_deepLTS_to_L2pyr * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supVIP    -> L2pyr
      do i = 1, num_supVIP_to_L2pyr
       j = map_supVIP_to_L2pyr(i,L) ! j = presynaptic cell
       k = com_supVIP_to_L2pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supVIP(j)  ! enumerate presyn. spikes
        presyntime = outtime_supVIP(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta)  ! 0.1 ms units, for otis

! GABA-A part
        dexparg = delta / tauGABA_supVIP_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L2pyr(k,L)  = gGABA_A_L2pyr(k,L) +
     &  gGABA_supVIP_to_L2pyr * z      
! end GABA-A part

c  k0 must be properly defined
      gGABA_B_L2pyr(k,L) = gGABA_B_L2pyr(k,L) +
     &   gGABAB_supVIP_to_L2pyr * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle placeholder5        -> L2pyr
      do i = 1, num_placeholder5_to_L2pyr
       j = map_placeholder5_to_L2pyr(i,L) ! j = presynaptic cell
       k = com_placeholder5_to_L2pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder5(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder5(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_placeholder5_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_L2pyr(k,L)  = gAMPA_L2pyr(k,L) +
     &  gAMPA_placeholder5_to_L2pyr * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_L2pyr(k,L) = gNMDA_L2pyr(k,L) +
     &  gNMDA_placeholder5_to_L2pyr * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_placeholder5_to_L2pyr
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_L2pyr(k,L) = gNMDA_L2pyr(k,L) +
     &  gNMDA_placeholder5_to_L2pyr * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_placeholder5_to_L2pyr
       if (gNMDA_L2pyr(k,L).gt.z)
     &  gNMDA_L2pyr(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle L3pyr  -> L2pyr
      do i = 1, num_L3pyr_to_L2pyr
       j = map_L3pyr_to_L2pyr(i,L) ! j = presynaptic cell
       k = com_L3pyr_to_L2pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L3pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L3pyr_to_L2pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_L2pyr(k,L)  = gAMPA_L2pyr(k,L) +
     &  gAMPA_L3pyr_to_L2pyr * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_L2pyr(k,L) = gNMDA_L2pyr(k,L) +
     &  gNMDA_L3pyr_to_L2pyr * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L3pyr_to_L2pyr
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_L2pyr(k,L) = gNMDA_L2pyr(k,L) +
     &  gNMDA_L3pyr_to_L2pyr * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L3pyr_to_L2pyr
       if (gNMDA_L2pyr(k,L).gt.z)
     &  gNMDA_L2pyr(k,L) = z
! end NMDA part

       end do ! m
      end do ! i

         end do
c End enumeration of L2pyr
       ENDIF ! if (mod(O,how_often).eq.0)...

! Define phasic currents to L2pyr cells, ectopic spikes,
! tonic synaptic conductances

      if (mod(O,200).eq.0) then
       call durand(seed,num_L2pyr,ranvec_L2pyr) 
!       do L = 1, num_L2pyr
        do L = firstcell, lastcell  ! Note
         if ((ranvec_L2pyr(L).gt.0.d0).and.
     &     (ranvec_L2pyr(L).le.noisepe_L2pyr)) then
c         curr_L2pyr(72,L) = 0.4d0
          curr_L2pyr(72,L) = 0.8d0
         else
          curr_L2pyr(72,L) = 0.d0
         endif 
        end do
      endif


! Call integration routine for L2pyr cells
       CALL INTEGRATE_L2pyr (O, time, num_L2pyr,
     &    V_L2pyr, curr_L2pyr,
     &    initialize, firstcell, lastcell,
     & gAMPA_L2pyr, gNMDA_L2pyr, gGABA_A_L2pyr,
     & gGABA_B_L2pyr, Mg, 
     & gapcon_L2pyr  ,totaxgj_L2pyr   ,gjtable_L2pyr, dt,
     &  chi_L2pyr,mnaf_L2pyr,mnap_L2pyr,
     &  hnaf_L2pyr,mkdr_L2pyr,mka_L2pyr,
     &  hka_L2pyr,mk2_L2pyr,hk2_L2pyr,
     &  mkm_L2pyr,mkc_L2pyr,mkahp_L2pyr,
     &  mcat_L2pyr,hcat_L2pyr,mcal_L2pyr,
     &  mar_L2pyr,field_sup ,field_deep, rel_axonshift_L2pyr)

  
       IF (mod(O,how_often).eq.0) then
! also field data                                     
c      do L = 1, num_L2pyr
       do L = firstcell, lastcell
        distal_axon_L2pyr (L-firstcell+1) = V_L2pyr (72,L)
       end do
  
           call mpi_allgather (distal_axon_L2pyr,
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_sup_local(1) = field_sup
        field_deep_local(1) = field_deep
           call mpi_allgather (field_sup_local,     
     &  1              , mpi_double_precision,
     &  field_sup_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_deep_local,     
     &  1              , mpi_double_precision,
     &  field_deep_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        ENDIF ! if (mod(O,how_often).eq.0) ....



! END thisno for L2pyr


c      ELSE IF (nodecell(thisno) .eq. 'placeholder1  ') THEN
       ELSE IF (nodecell(thisno) .eq. 'supintern   ') THEN
c placeholder1

c Determine which particular cells this node will be concerned with.
c         i = place (thisno)
          firstcell = 1 
          lastcell =  num_placeholder1 

        IF (mod(O,how_often).eq.0) then
c 1st set placeholder1 synaptic conductances to 0:

          do i = 1, numcomp_placeholder1
          do j = firstcell, lastcell
         gAMPA_placeholder1(i,j)     = 0.d0
         gNMDA_placeholder1(i,j)     = 0.d0
         gGABA_A_placeholder1(i,j)   = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle L2pyr   -> placeholder1
      do i = 1, num_L2pyr_to_placeholder1  
       j = map_L2pyr_to_placeholder1(i,L) ! j = presynaptic cell
       k = com_L2pyr_to_placeholder1(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L2pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L2pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L2pyr_to_placeholder1  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder1(k,L)  = gAMPA_placeholder1(k,L) +
     &  gAMPA_L2pyr_to_placeholder1 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder1(k,L) = gNMDA_placeholder1(k,L) +
     &  gNMDA_L2pyr_to_placeholder1 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L2pyr_to_placeholder1  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder1(k,L) = gNMDA_placeholder1(k,L) +
     &  gNMDA_L2pyr_to_placeholder1 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L2pyr_to_placeholder1  
       if (gNMDA_placeholder1(k,L).gt.z)
     &  gNMDA_placeholder1(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle placeholder1    -> placeholder1
      do i = 1, num_placeholder1_to_placeholder1  
       j = map_placeholder1_to_placeholder1(i,L) ! j = presynaptic cell
       k = com_placeholder1_to_placeholder1(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder1(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder1(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder1_to_placeholder1  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_placeholder1(k,L)  = gGABA_A_placeholder1(k,L) +
     &  gGABA_placeholder1_to_placeholder1 * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle supng      -> placeholder1 
      do i = 1, num_supng_to_placeholder1 
       j = map_supng_to_placeholder1 (i,L) ! j = presynaptic cell
       k = com_supng_to_placeholder1 (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supng(j)  ! enumerate presyn. spikes
        presyntime = outtime_supng(m,j)
        delta = time - presyntime

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_supng_to_placeholder1 
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_placeholder1 (k,L)  = gGABA_A_placeholder1 (k,L) +
     &  gGABA_supng_to_placeholder1  * z      
! end GABA-A part

       end do ! m
      end do ! i



c Handle placeholder3     -> placeholder1
      do i = 1, num_placeholder3_to_placeholder1  
       j = map_placeholder3_to_placeholder1(i,L) ! j = presynaptic cell
       k = com_placeholder3_to_placeholder1(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder3(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder3(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder3_to_placeholder1  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_placeholder1(k,L)  = gGABA_A_placeholder1(k,L) +
     &  gGABA_placeholder3_to_placeholder1 * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle LOT  -> placeholder1
      do i = 1, num_LOT_to_placeholder1  
       j = map_LOT_to_placeholder1(i,L) ! j = presynaptic cell
       k = com_LOT_to_placeholder1(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_LOT(j)  ! enumerate presyn. spikes
        presyntime = outtime_LOT(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_LOT_to_placeholder1  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder1(k,L)  = gAMPA_placeholder1(k,L) +
     &  gAMPA_LOT_to_placeholder1 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder1(k,L) = gNMDA_placeholder1(k,L) +
     &  gNMDA_LOT_to_placeholder1 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_LOT_to_placeholder1  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder1(k,L) = gNMDA_placeholder1(k,L) +
     &  gNMDA_LOT_to_placeholder1 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_LOT_to_placeholder1  
       if (gNMDA_placeholder1(k,L).gt.z)
     &  gNMDA_placeholder1(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle semilunar     -> placeholder1
      do i = 1, num_semilunar_to_placeholder1  
       j = map_semilunar_to_placeholder1(i,L) ! j = presynaptic cell
       k = com_semilunar_to_placeholder1(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_semilunar(j)  ! enumerate presyn. spikes
        presyntime = outtime_semilunar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_semilunar_to_placeholder1  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder1(k,L)  = gAMPA_placeholder1(k,L) +
     &  gAMPA_semilunar_to_placeholder1 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder1(k,L) = gNMDA_placeholder1(k,L) +
     &  gNMDA_semilunar_to_placeholder1 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_semilunar_to_placeholder1  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder1(k,L) = gNMDA_placeholder1(k,L) +
     &  gNMDA_semilunar_to_placeholder1 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_semilunar_to_placeholder1  
       if (gNMDA_placeholder1(k,L).gt.z)
     &  gNMDA_placeholder1(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle multipolar     -> placeholder1
      do i = 1, num_multipolar_to_placeholder1  
       j = map_multipolar_to_placeholder1(i,L) ! j = presynaptic cell
       k = com_multipolar_to_placeholder1(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_multipolar(j)  ! enumerate presyn. spikes
        presyntime = outtime_multipolar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_multipolar_to_placeholder1  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder1(k,L)  = gAMPA_placeholder1(k,L) +
     &  gAMPA_multipolar_to_placeholder1 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder1(k,L) = gNMDA_placeholder1(k,L) +
     &  gNMDA_multipolar_to_placeholder1 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_multipolar_to_placeholder1  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder1(k,L) = gNMDA_placeholder1(k,L) +
     &  gNMDA_multipolar_to_placeholder1 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_multipolar_to_placeholder1  
       if (gNMDA_placeholder1(k,L).gt.z)
     &  gNMDA_placeholder1(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supVIP    -> placeholder1
      do i = 1, num_supVIP_to_placeholder1  
       j = map_supVIP_to_placeholder1(i,L) ! j = presynaptic cell
       k = com_supVIP_to_placeholder1(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supVIP(j)  ! enumerate presyn. spikes
        presyntime = outtime_supVIP(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supVIP_to_placeholder1  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_placeholder1(k,L)  = gGABA_A_placeholder1(k,L) +
     &  gGABA_supVIP_to_placeholder1 * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle placeholder5    -> placeholder1
      do i = 1, num_placeholder5_to_placeholder1  
       j = map_placeholder5_to_placeholder1(i,L) ! j = presynaptic cell
       k = com_placeholder5_to_placeholder1(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder5(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder5(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_placeholder5_to_placeholder1  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder1(k,L)  = gAMPA_placeholder1(k,L) +
     &  gAMPA_placeholder5_to_placeholder1 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder1(k,L) = gNMDA_placeholder1(k,L) +
     &  gNMDA_placeholder5_to_placeholder1 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_placeholder5_to_placeholder1  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder1(k,L) = gNMDA_placeholder1(k,L) +
     &  gNMDA_placeholder5_to_placeholder1 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_placeholder5_to_placeholder1  
       if (gNMDA_placeholder1(k,L).gt.z)
     &  gNMDA_placeholder1(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle L3pyr  -> placeholder1
      do i = 1, num_L3pyr_to_placeholder1  
       j = map_L3pyr_to_placeholder1(i,L) ! j = presynaptic cell
       k = com_L3pyr_to_placeholder1(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L3pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L3pyr_to_placeholder1  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder1(k,L)  = gAMPA_placeholder1(k,L) +
     &  gAMPA_L3pyr_to_placeholder1 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder1(k,L) = gNMDA_placeholder1(k,L) +
     &  gNMDA_L3pyr_to_placeholder1 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L3pyr_to_placeholder1  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder1(k,L) = gNMDA_placeholder1(k,L) +
     &  gNMDA_L3pyr_to_placeholder1 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L3pyr_to_placeholder1  
       if (gNMDA_placeholder1(k,L).gt.z)
     &  gNMDA_placeholder1(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of placeholder1  
        ENDIF  ! if (mod(O,how_often).eq.0) ....

! Define currents to placeholder1   cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for placeholder1   cells
       CALL INTEGRATE_placeholder1 (O, time, num_placeholder1 ,
     &    V_placeholder1 , curr_placeholder1 ,
     $    initialize, firstcell, lastcell,
     & gAMPA_placeholder1 , gNMDA_placeholder1 , gGABA_A_placeholder1 ,
     & Mg, 
     & gapcon_placeholder1,totSDgj_placeholder1,gjtable_placeholder1,dt,
     &  chi_placeholder1,mnaf_placeholder1,mnap_placeholder1,
     &  hnaf_placeholder1,mkdr_placeholder1,mka_placeholder1,
     &  hka_placeholder1,mk2_placeholder1,hk2_placeholder1,
     &  mkm_placeholder1,mkc_placeholder1,mkahp_placeholder1,
     &  mcat_placeholder1,hcat_placeholder1,mcal_placeholder1,
     &  mar_placeholder1)


      IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_placeholder1  
       do L = firstcell, lastcell
c       distal_axon_placeholder1   (L-firstcell+1) = V_placeholder1   (59,L)
        distal_axon_supintern (L-firstcell+1) = V_placeholder1   (59,L)
       end do
  
c          call mpi_allgather (distal_axon_placeholder1,
c    &  maxcellspernode, mpi_double_precision,
c    &  distal_axon_global,maxcellspernode,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)

        field_sup_local(1) = 0.d0     
        field_deep_local(1) = 0.d0     
c          call mpi_allgather (field_sup_local,     
c    &  1              , mpi_double_precision,
c    &  field_sup_global  , 1             ,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)
c          call mpi_allgather (field_deep_local,     
c    &  1              , mpi_double_precision,
c    &  field_deep_global  , 1             ,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)
  
           ENDIF  ! if (mod(O,how_often).eq.0) ....

! END thisno for placeholder1

c      ELSE IF (nodecell(thisno) .eq. 'supng    ') THEN
c supng  

c Determine which particular cells this node will be concerned with.
c         i = place (thisno)
          firstcell = 1 
          lastcell =  num_supng 

        IF (mod(O,how_often).eq.0) then
c 1st set supng   synaptic conductances to 0:

          do i = 1, numcomp_placeholder1
          do j = firstcell, lastcell
         gAMPA_supng  (i,j)     = 0.d0
         gNMDA_supng  (i,j)     = 0.d0
         gGABA_A_supng  (i,j)   = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle L2pyr   -> supng  
      do i = 1, num_L2pyr_to_supng    
       j = map_L2pyr_to_supng  (i,L) ! j = presynaptic cell
       k = com_L2pyr_to_supng  (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L2pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L2pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L2pyr_to_supng    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supng  (k,L)  = gAMPA_supng  (k,L) +
     &  gAMPA_L2pyr_to_supng   * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supng  (k,L) = gNMDA_supng  (k,L) +
     &  gNMDA_L2pyr_to_supng   * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L2pyr_to_supng    
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supng  (k,L) = gNMDA_supng  (k,L) +
     &  gNMDA_L2pyr_to_supng   * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L2pyr_to_supng    
       if (gNMDA_supng  (k,L).gt.z)
     &  gNMDA_supng  (k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle placeholder1    -> supng  
      do i = 1, num_placeholder1_to_supng    
       j = map_placeholder1_to_supng  (i,L) ! j = presynaptic cell
       k = com_placeholder1_to_supng  (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder1(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder1(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder1_to_supng    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supng  (k,L)  = gGABA_A_supng  (k,L) +
     &  gGABA_placeholder1_to_supng   * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle supng      -> supng   
      do i = 1, num_supng_to_supng   
       j = map_supng_to_supng   (i,L) ! j = presynaptic cell
       k = com_supng_to_supng   (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supng(j)  ! enumerate presyn. spikes
        presyntime = outtime_supng(m,j)
        delta = time - presyntime

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_supng_to_supng   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supng   (k,L)  = gGABA_A_supng   (k,L) +
     &  gGABA_supng_to_supng    * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle supVIP     -> supng   
      do i = 1, num_supVIP_to_supng   
       j = map_supVIP_to_supng   (i,L) ! j = presynaptic cell
       k = com_supVIP_to_supng   (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supVIP(j)  ! enumerate presyn. spikes
        presyntime = outtime_supVIP(m,j)
        delta = time - presyntime

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_supVIP_to_supng   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supng   (k,L)  = gGABA_A_supng   (k,L) +
     &  gGABA_supVIP_to_supng    * z      
! end GABA-A part

       end do ! m
      end do ! i



c Handle LOT  -> supng 
      do i = 1, num_LOT_to_supng     
       j = map_LOT_to_supng (i,L) ! j = presynaptic cell
       k = com_LOT_to_supng (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_LOT(j)  ! enumerate presyn. spikes
        presyntime = outtime_LOT(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_LOT_to_supng    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supng (k,L)  = gAMPA_supng (k,L) +
     &  gAMPA_LOT_to_supng  * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supng (k,L) = gNMDA_supng (k,L) +
     &  gNMDA_LOT_to_supng  * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_LOT_to_supng    
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supng (k,L) = gNMDA_supng (k,L) +
     &  gNMDA_LOT_to_supng  * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_LOT_to_supng   
       if (gNMDA_supng (k,L).gt.z)
     &  gNMDA_supng (k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle placeholder5    -> supng  
      do i = 1, num_placeholder5_to_supng    
       j = map_placeholder5_to_supng  (i,L) ! j = presynaptic cell
       k = com_placeholder5_to_supng  (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder5(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder5(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_placeholder5_to_supng    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supng  (k,L)  = gAMPA_supng  (k,L) +
     &  gAMPA_placeholder5_to_supng   * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supng  (k,L) = gNMDA_supng  (k,L) +
     &  gNMDA_placeholder5_to_supng   * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_placeholder5_to_supng    
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supng  (k,L) = gNMDA_supng  (k,L) +
     &  gNMDA_placeholder5_to_supng   * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_placeholder5_to_supng    
       if (gNMDA_supng  (k,L).gt.z)
     &  gNMDA_supng  (k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


         end do
c End enumeration of supng    
        ENDIF  ! if (mod(O,how_often).eq.0) ....

! Define currents to supng     cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for supng     cells
       CALL INTEGRATE_supng    (O, time, num_supng   ,
     &    V_supng   , curr_supng   ,
     $    initialize, firstcell, lastcell,
     & gAMPA_supng   , gNMDA_supng   , gGABA_A_supng   ,
     & Mg, 
     & gapcon_supng     ,totSDgj_supng      ,gjtable_supng   , dt,
     &  chi_supng  ,mnaf_supng  ,mnap_supng  ,
     &  hnaf_supng  ,mkdr_supng  ,mka_supng  ,
     &  hka_supng  ,mk2_supng  ,hk2_supng  ,
     &  mkm_supng  ,mkc_supng  ,mkahp_supng  ,
     &  mcat_supng  ,hcat_supng  ,mcal_supng  ,
     &  mar_supng  )


      IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_placeholder1  
       do L = firstcell, lastcell
        distal_axon_supintern  (L + 400) =  V_supng     (59,L)
       end do
  
c          call mpi_allgather (distal_axon_supng  ,
c    &  maxcellspernode, mpi_double_precision,
c    &  distal_axon_global,maxcellspernode,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)

        field_sup_local(1) = 0.d0     
        field_deep_local(1) = 0.d0     
c          call mpi_allgather (field_sup_local,     
c    &  1              , mpi_double_precision,
c    &  field_sup_global  , 1             ,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)
c          call mpi_allgather (field_deep_local,     
c    &  1              , mpi_double_precision,
c    &  field_deep_global  , 1             ,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)
  
           ENDIF  ! if (mod(O,how_often).eq.0) ....

! END thisno for supng  

c      ELSE IF (THISNO.EQ.3) THEN
c      ELSE IF (nodecell(thisno) .eq. 'placeholder2  ') THEN
c placeholder2

c Determine which particular cells this node will be concerned with.
c         i = place (thisno)
          firstcell = 1 
          lastcell =  num_placeholder2 

         IF (mod(O,how_often).eq.0) then
c 1st set placeholder2 synaptic conductances to 0:

          do i = 1, numcomp_placeholder2
          do j = firstcell, lastcell
         gAMPA_placeholder2(i,j)     = 0.d0
         gNMDA_placeholder2(i,j)     = 0.d0
         gGABA_A_placeholder2(i,j)   = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle L2pyr   -> placeholder2
      do i = 1, num_L2pyr_to_placeholder2  
       j = map_L2pyr_to_placeholder2(i,L) ! j = presynaptic cell
       k = com_L2pyr_to_placeholder2(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L2pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L2pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L2pyr_to_placeholder2  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder2(k,L)  = gAMPA_placeholder2(k,L) +
     &  gAMPA_L2pyr_to_placeholder2 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder2(k,L) = gNMDA_placeholder2(k,L) +
     &  gNMDA_L2pyr_to_placeholder2 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L2pyr_to_placeholder2  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder2(k,L) = gNMDA_placeholder2(k,L) +
     &  gNMDA_L2pyr_to_placeholder2 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L2pyr_to_placeholder2  
       if (gNMDA_placeholder2(k,L).gt.z)
     &  gNMDA_placeholder2(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle placeholder1    -> placeholder2
      do i = 1, num_placeholder1_to_placeholder2  
       j = map_placeholder1_to_placeholder2(i,L) ! j = presynaptic cell
       k = com_placeholder1_to_placeholder2(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder1(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder1(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder1_to_placeholder2  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_placeholder2(k,L)  = gGABA_A_placeholder2(k,L) +
     &  gGABA_placeholder1_to_placeholder2 * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle placeholder3     -> placeholder2
      do i = 1, num_placeholder3_to_placeholder2  
       j = map_placeholder3_to_placeholder2(i,L) ! j = presynaptic cell
       k = com_placeholder3_to_placeholder2(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder3(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder3(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder3_to_placeholder2  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_placeholder2(k,L)  = gGABA_A_placeholder2(k,L) +
     &  gGABA_placeholder3_to_placeholder2 * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle LOT  -> placeholder2
      do i = 1, num_LOT_to_placeholder2  
       j = map_LOT_to_placeholder2(i,L) ! j = presynaptic cell
       k = com_LOT_to_placeholder2(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_LOT(j)  ! enumerate presyn. spikes
        presyntime = outtime_LOT(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_LOT_to_placeholder2  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder2(k,L)  = gAMPA_placeholder2(k,L) +
     &  gAMPA_LOT_to_placeholder2 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder2(k,L) = gNMDA_placeholder2(k,L) +
     &  gNMDA_LOT_to_placeholder2 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_LOT_to_placeholder2  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder2(k,L) = gNMDA_placeholder2(k,L) +
     &  gNMDA_LOT_to_placeholder2 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_LOT_to_placeholder2  
       if (gNMDA_placeholder2(k,L).gt.z)
     &  gNMDA_placeholder2(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle semilunar     -> placeholder2
      do i = 1, num_semilunar_to_placeholder2  
       j = map_semilunar_to_placeholder2(i,L) ! j = presynaptic cell
       k = com_semilunar_to_placeholder2(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_semilunar(j)  ! enumerate presyn. spikes
        presyntime = outtime_semilunar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_semilunar_to_placeholder2  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder2(k,L)  = gAMPA_placeholder2(k,L) +
     &  gAMPA_semilunar_to_placeholder2 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder2(k,L) = gNMDA_placeholder2(k,L) +
     &  gNMDA_semilunar_to_placeholder2 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_semilunar_to_placeholder2  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder2(k,L) = gNMDA_placeholder2(k,L) +
     &  gNMDA_semilunar_to_placeholder2 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_semilunar_to_placeholder2  
       if (gNMDA_placeholder2(k,L).gt.z)
     &  gNMDA_placeholder2(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle multipolar     -> placeholder2
      do i = 1, num_multipolar_to_placeholder2  
       j = map_multipolar_to_placeholder2(i,L) ! j = presynaptic cell
       k = com_multipolar_to_placeholder2(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_multipolar(j)  ! enumerate presyn. spikes
        presyntime = outtime_multipolar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_multipolar_to_placeholder2  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder2(k,L)  = gAMPA_placeholder2(k,L) +
     &  gAMPA_multipolar_to_placeholder2 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder2(k,L) = gNMDA_placeholder2(k,L) +
     &  gNMDA_multipolar_to_placeholder2 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_multipolar_to_placeholder2  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder2(k,L) = gNMDA_placeholder2(k,L) +
     &  gNMDA_multipolar_to_placeholder2 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_multipolar_to_placeholder2  
       if (gNMDA_placeholder2(k,L).gt.z)
     &  gNMDA_placeholder2(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supVIP    -> placeholder2
      do i = 1, num_supVIP_to_placeholder2  
       j = map_supVIP_to_placeholder2(i,L) ! j = presynaptic cell
       k = com_supVIP_to_placeholder2(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supVIP(j)  ! enumerate presyn. spikes
        presyntime = outtime_supVIP(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supVIP_to_placeholder2  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_placeholder2(k,L)  = gGABA_A_placeholder2(k,L) +
     &  gGABA_supVIP_to_placeholder2 * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle placeholder5        -> placeholder2
      do i = 1, num_placeholder5_to_placeholder2  
       j = map_placeholder5_to_placeholder2(i,L) ! j = presynaptic cell
       k = com_placeholder5_to_placeholder2(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder5(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder5(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_placeholder5_to_placeholder2  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder2(k,L)  = gAMPA_placeholder2(k,L) +
     &  gAMPA_placeholder5_to_placeholder2 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder2(k,L) = gNMDA_placeholder2(k,L) +
     &  gNMDA_placeholder5_to_placeholder2 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_placeholder5_to_placeholder2  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder2(k,L) = gNMDA_placeholder2(k,L) +
     &  gNMDA_placeholder5_to_placeholder2 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_placeholder5_to_placeholder2  
       if (gNMDA_placeholder2(k,L).gt.z)
     &  gNMDA_placeholder2(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle L3pyr  -> placeholder2
      do i = 1, num_L3pyr_to_placeholder2  
       j = map_L3pyr_to_placeholder2(i,L) ! j = presynaptic cell
       k = com_L3pyr_to_placeholder2(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L3pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L3pyr_to_placeholder2  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder2(k,L)  = gAMPA_placeholder2(k,L) +
     &  gAMPA_L3pyr_to_placeholder2 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder2(k,L) = gNMDA_placeholder2(k,L) +
     &  gNMDA_L3pyr_to_placeholder2 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L3pyr_to_placeholder2  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder2(k,L) = gNMDA_placeholder2(k,L) +
     &  gNMDA_L3pyr_to_placeholder2 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L3pyr_to_placeholder2  
       if (gNMDA_placeholder2(k,L).gt.z)
     &  gNMDA_placeholder2(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of placeholder2  
         ENDIF  ! if (mod(O,how_often).eq.0) ...


! Define currents to placeholder2   cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for placeholder2   cells
       CALL INTEGRATE_placeholder2 (O, time, num_placeholder2 ,
     &    V_placeholder2 , curr_placeholder2 ,
     &    initialize, firstcell, lastcell,
     & gAMPA_placeholder2 , gNMDA_placeholder2 , gGABA_A_placeholder2 ,
     & Mg, 
     & gapcon_placeholder2,totSDgj_placeholder2,gjtable_placeholder2,dt,
     &  chi_placeholder2,mnaf_placeholder2,mnap_placeholder2,
     &  hnaf_placeholder2,mkdr_placeholder2,mka_placeholder2,
     &  hka_placeholder2,mk2_placeholder2,hk2_placeholder2,
     &  mkm_placeholder2,mkc_placeholder2,mkahp_placeholder2,
     &  mcat_placeholder2,hcat_placeholder2,mcal_placeholder2,
     &  mar_placeholder2)


        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_placeholder2  
       do L = firstcell, lastcell
        distal_axon_supintern (L + 100     ) = V_placeholder2   (59,L)
       end do
  
c          call mpi_allgather (distal_axon_placeholder2, 
c    &  maxcellspernode, mpi_double_precision,
c    &  distal_axon_global,maxcellspernode,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)

        field_sup_local(1) = 0.d0     
        field_deep_local(1) = 0.d0     
c          call mpi_allgather (field_sup_local,     
c    &  1              , mpi_double_precision,
c    &  field_sup_global  , 1             ,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)
c          call mpi_allgather (field_deep_local,     
c    &  1              , mpi_double_precision,
c    &  field_deep_global  , 1             ,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)
  
             ENDIF !  if (mod(O,how_often).eq.0) ...

! END thisno for placeholder2

c      ELSE IF (THISNO.EQ.4) THEN
c      ELSE IF (nodecell(thisno) .eq. 'placeholder3   ') THEN
c placeholder3

c Determine which particular cells this node will be concerned with.
c         i = place (thisno)
          firstcell = 1 
          lastcell = num_placeholder3                            

          IF (mod(O,how_often).eq.0) then
c 1st set placeholder3  synaptic conductances to 0:

          do i = 1, numcomp_placeholder3
          do j = firstcell, lastcell
         gAMPA_placeholder3(i,j)      = 0.d0
         gNMDA_placeholder3(i,j)      = 0.d0
         gGABA_A_placeholder3(i,j)    = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle L2pyr   -> placeholder3
      do i = 1, num_L2pyr_to_placeholder3   
       j = map_L2pyr_to_placeholder3(i,L) ! j = presynaptic cell
       k = com_L2pyr_to_placeholder3(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L2pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L2pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L2pyr_to_placeholder3  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder3(k,L)  = gAMPA_placeholder3(k,L) +
     &  gAMPA_L2pyr_to_placeholder3 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder3(k,L) = gNMDA_placeholder3(k,L) +
     &  gNMDA_L2pyr_to_placeholder3 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L2pyr_to_placeholder3  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder3(k,L) = gNMDA_placeholder3(k,L) +
     &  gNMDA_L2pyr_to_placeholder3 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L2pyr_to_placeholder3  
       if (gNMDA_placeholder3(k,L).gt.z)
     &  gNMDA_placeholder3(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle placeholder1    -> placeholder3
      do i = 1, num_placeholder1_to_placeholder3  
       j = map_placeholder1_to_placeholder3(i,L) ! j = presynaptic cell
       k = com_placeholder1_to_placeholder3(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder1(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder1(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder1_to_placeholder3  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_placeholder3(k,L)  = gGABA_A_placeholder3(k,L) +
     &  gGABA_placeholder1_to_placeholder3 * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle placeholder3     -> placeholder3
      do i = 1, num_placeholder3_to_placeholder3  
       j = map_placeholder3_to_placeholder3(i,L) ! j = presynaptic cell
       k = com_placeholder3_to_placeholder3(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder3(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder3(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder3_to_placeholder3  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_placeholder3(k,L)  = gGABA_A_placeholder3(k,L) +
     &  gGABA_placeholder3_to_placeholder3 * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle LOT  -> placeholder3
      do i = 1, num_LOT_to_placeholder3  
       j = map_LOT_to_placeholder3(i,L) ! j = presynaptic cell
       k = com_LOT_to_placeholder3(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_LOT(j)  ! enumerate presyn. spikes
        presyntime = outtime_LOT(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_LOT_to_placeholder3  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder3(k,L)  = gAMPA_placeholder3(k,L) +
     &  gAMPA_LOT_to_placeholder3 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder3(k,L) = gNMDA_placeholder3(k,L) +
     &  gNMDA_LOT_to_placeholder3 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_LOT_to_placeholder3  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder3(k,L) = gNMDA_placeholder3(k,L) +
     &  gNMDA_LOT_to_placeholder3 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_LOT_to_placeholder3  
       if (gNMDA_placeholder3(k,L).gt.z)
     &  gNMDA_placeholder3(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle semilunar     -> placeholder3
      do i = 1, num_semilunar_to_placeholder3  
       j = map_semilunar_to_placeholder3(i,L) ! j = presynaptic cell
       k = com_semilunar_to_placeholder3(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_semilunar(j)  ! enumerate presyn. spikes
        presyntime = outtime_semilunar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_semilunar_to_placeholder3  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder3(k,L)  = gAMPA_placeholder3(k,L) +
     &  gAMPA_semilunar_to_placeholder3 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder3(k,L) = gNMDA_placeholder3(k,L) +
     &  gNMDA_semilunar_to_placeholder3 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_semilunar_to_placeholder3  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder3(k,L) = gNMDA_placeholder3(k,L) +
     &  gNMDA_semilunar_to_placeholder3 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_semilunar_to_placeholder3  
       if (gNMDA_placeholder3(k,L).gt.z)
     &  gNMDA_placeholder3(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle multipolar     -> placeholder3
      do i = 1, num_multipolar_to_placeholder3  
       j = map_multipolar_to_placeholder3(i,L) ! j = presynaptic cell
       k = com_multipolar_to_placeholder3(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_multipolar(j)  ! enumerate presyn. spikes
        presyntime = outtime_multipolar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_multipolar_to_placeholder3  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder3(k,L)  = gAMPA_placeholder3(k,L) +
     &  gAMPA_multipolar_to_placeholder3 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder3(k,L) = gNMDA_placeholder3(k,L) +
     &  gNMDA_multipolar_to_placeholder3 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_multipolar_to_placeholder3  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder3(k,L) = gNMDA_placeholder3(k,L) +
     &  gNMDA_multipolar_to_placeholder3 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_multipolar_to_placeholder3  
       if (gNMDA_placeholder3(k,L).gt.z)
     &  gNMDA_placeholder3(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supVIP    -> placeholder3
      do i = 1, num_supVIP_to_placeholder3   
       j = map_supVIP_to_placeholder3(i,L) ! j = presynaptic cell
       k = com_supVIP_to_placeholder3(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supVIP(j)  ! enumerate presyn. spikes
        presyntime = outtime_supVIP(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supVIP_to_placeholder3   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_placeholder3(k,L)  = gGABA_A_placeholder3(k,L) +
     &  gGABA_supVIP_to_placeholder3 * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle L3pyr  -> placeholder3
      do i = 1, num_L3pyr_to_placeholder3  
       j = map_L3pyr_to_placeholder3(i,L) ! j = presynaptic cell
       k = com_L3pyr_to_placeholder3(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L3pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L3pyr_to_placeholder3  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder3(k,L)  = gAMPA_placeholder3(k,L) +
     &  gAMPA_L3pyr_to_placeholder3 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder3(k,L) = gNMDA_placeholder3(k,L) +
     &  gNMDA_L3pyr_to_placeholder3 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L3pyr_to_placeholder3  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder3(k,L) = gNMDA_placeholder3(k,L) +
     &  gNMDA_L3pyr_to_placeholder3 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L3pyr_to_placeholder3  
       if (gNMDA_placeholder3(k,L).gt.z)
     &  gNMDA_placeholder3(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of placeholder3   
        ENDIF  ! if (mod(O,how_often).eq.0) ...

! Define currents to placeholder3    cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for placeholder3    cells

       CALL INTEGRATE_placeholder3  (O, time, num_placeholder3  ,
     &    V_placeholder3  , curr_placeholder3  ,
     &    initialize, firstcell, lastcell,
     & gAMPA_placeholder3,gNMDA_placeholder3  , gGABA_A_placeholder3  ,
     & Mg, 
     & gapcon_placeholder3,totSDgj_placeholder3,gjtable_placeholder3,dt,
     &  chi_placeholder3,mnaf_placeholder3,mnap_placeholder3,
     &  hnaf_placeholder3,mkdr_placeholder3,mka_placeholder3,
     &  hka_placeholder3,mk2_placeholder3,hk2_placeholder3,
     &  mkm_placeholder3,mkc_placeholder3,mkahp_placeholder3,
     &  mcat_placeholder3,hcat_placeholder3,mcal_placeholder3,
     &  mar_placeholder3)


        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
       do L = 1, num_placeholder3   
        distal_axon_supintern (L + 200    ) = V_placeholder3    (59,L)
       end do
  
c          call mpi_allgather (distal_axon_placeholder3,   
c    &  maxcellspernode, mpi_double_precision,
c    &  distal_axon_global,maxcellspernode,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)

        field_sup_local(1) = 0.d0     
        field_deep_local(1) = 0.d0     
c          call mpi_allgather (field_sup_local,     
c    &  1              , mpi_double_precision,
c    &  field_sup_global  , 1             ,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)
c          call mpi_allgather (field_deep_local,     
c    &  1              , mpi_double_precision,
c    &  field_deep_global  , 1             ,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)
  
         ENDIF  ! if (mod(O,how_often).eq.0) ...

! END thisno for placeholder3

c      ELSE IF (nodecell(thisno) .eq. 'supVIP  ') THEN
c supVIP

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 
          lastcell = num_supVIP                            

       IF (mod(O,how_often).eq.0) then
c 1st set supVIP   synaptic conductances to 0:

          do i = 1, numcomp_supVIP
          do j = firstcell, lastcell
         gAMPA_supVIP(i,j)     = 0.d0
         gNMDA_supVIP(i,j)     = 0.d0
         gGABA_A_supVIP(i,j)   = 0.d0 
          end do
          end do

         do L = firstcell, lastcell
c Handle L2pyr   -> supVIP
      do i = 1, num_L2pyr_to_supVIP   
       j = map_L2pyr_to_supVIP(i,L) ! j = presynaptic cell
       k = com_L2pyr_to_supVIP(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L2pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L2pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L2pyr_to_supVIP  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supVIP(k,L)  = gAMPA_supVIP(k,L) +
     &  gAMPA_L2pyr_to_supVIP * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supVIP(k,L) = gNMDA_supVIP(k,L) +
     &  gNMDA_L2pyr_to_supVIP * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L2pyr_to_supVIP  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supVIP(k,L) = gNMDA_supVIP(k,L) +
     &  gNMDA_L2pyr_to_supVIP * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L2pyr_to_supVIP  
       if (gNMDA_supVIP(k,L).gt.z)
     &  gNMDA_supVIP(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle placeholder3     -> supVIP
      do i = 1, num_placeholder3_to_supVIP     
       j = map_placeholder3_to_supVIP(i,L) ! j = presynaptic cell
       k = com_placeholder3_to_supVIP(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder3(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder3(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder3_to_supVIP     
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supVIP(k,L)  = gGABA_A_supVIP(k,L) +
     &  gGABA_placeholder3_to_supVIP * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle LOT  -> supVIP
      do i = 1, num_LOT_to_supVIP    
       j = map_LOT_to_supVIP(i,L) ! j = presynaptic cell
       k = com_LOT_to_supVIP(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_LOT(j)  ! enumerate presyn. spikes
        presyntime = outtime_LOT(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_LOT_to_supVIP   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supVIP(k,L)  = gAMPA_supVIP(k,L) +
     &  gAMPA_LOT_to_supVIP * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supVIP(k,L) = gNMDA_supVIP(k,L) +
     &  gNMDA_LOT_to_supVIP * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_LOT_to_supVIP   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supVIP(k,L) = gNMDA_supVIP(k,L) +
     &  gNMDA_LOT_to_supVIP * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_LOT_to_supVIP  
       if (gNMDA_supVIP(k,L).gt.z)
     &  gNMDA_supVIP(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle semilunar     -> supVIP
      do i = 1, num_semilunar_to_supVIP    
       j = map_semilunar_to_supVIP(i,L) ! j = presynaptic cell
       k = com_semilunar_to_supVIP(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_semilunar(j)  ! enumerate presyn. spikes
        presyntime = outtime_semilunar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_semilunar_to_supVIP   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supVIP(k,L)  = gAMPA_supVIP(k,L) +
     &  gAMPA_semilunar_to_supVIP * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supVIP(k,L) = gNMDA_supVIP(k,L) +
     &  gNMDA_semilunar_to_supVIP * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_semilunar_to_supVIP   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supVIP(k,L) = gNMDA_supVIP(k,L) +
     &  gNMDA_semilunar_to_supVIP * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_semilunar_to_supVIP   
       if (gNMDA_supVIP(k,L).gt.z)
     &  gNMDA_supVIP(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle multipolar     -> supVIP
      do i = 1, num_multipolar_to_supVIP    
       j = map_multipolar_to_supVIP(i,L) ! j = presynaptic cell
       k = com_multipolar_to_supVIP(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_multipolar(j)  ! enumerate presyn. spikes
        presyntime = outtime_multipolar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_multipolar_to_supVIP   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supVIP(k,L)  = gAMPA_supVIP(k,L) +
     &  gAMPA_multipolar_to_supVIP * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supVIP(k,L) = gNMDA_supVIP(k,L) +
     &  gNMDA_multipolar_to_supVIP * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_multipolar_to_supVIP   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supVIP(k,L) = gNMDA_supVIP(k,L) +
     &  gNMDA_multipolar_to_supVIP * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_multipolar_to_supVIP   
       if (gNMDA_supVIP(k,L).gt.z)
     &  gNMDA_supVIP(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask   -> supVIP
      do i = 1, num_deepbask_to_supVIP     
       j = map_deepbask_to_supVIP(i,L) ! j = presynaptic cell
       k = com_deepbask_to_supVIP(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_supVIP     
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supVIP(k,L)  = gGABA_A_supVIP(k,L) +
     &  gGABA_deepbask_to_supVIP * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supVIP    -> supVIP
      do i = 1, num_supVIP_to_supVIP     
       j = map_supVIP_to_supVIP(i,L) ! j = presynaptic cell
       k = com_supVIP_to_supVIP(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supVIP(j)  ! enumerate presyn. spikes
        presyntime = outtime_supVIP(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supVIP_to_supVIP     
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_supVIP(k,L)  = gGABA_A_supVIP(k,L) +
     &  gGABA_supVIP_to_supVIP * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle L3pyr  -> supVIP
      do i = 1, num_L3pyr_to_supVIP
       j = map_L3pyr_to_supVIP(i,L) ! j = presynaptic cell
       k = com_L3pyr_to_supVIP(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L3pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L3pyr_to_supVIP
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_supVIP(k,L)  = gAMPA_supVIP(k,L) +
     &  gAMPA_L3pyr_to_supVIP * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_supVIP(k,L) = gNMDA_supVIP(k,L) +
     &  gNMDA_L3pyr_to_supVIP * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L3pyr_to_supVIP 
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_supVIP(k,L) = gNMDA_supVIP(k,L) +
     &  gNMDA_L3pyr_to_supVIP * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L3pyr_to_supVIP
       if (gNMDA_supVIP(k,L).gt.z)
     &  gNMDA_supVIP(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of supVIP     
         ENDIF  !  if (mod(O,how_often).eq.0) ...

! Define currents to supVIP      cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for supVIP      cells
       CALL INTEGRATE_supVIP  (O, time, num_supVIP  ,
     &    V_supVIP  , curr_supVIP  ,
     & initialize, firstcell, lastcell,
     & gAMPA_supVIP  , gNMDA_supVIP  , gGABA_A_supVIP  ,
     & Mg, 
     & gapcon_supVIP  ,totSDgj_supVIP  ,gjtable_supVIP  , dt,
     &  chi_supVIP,mnaf_supVIP,mnap_supVIP,
     &  hnaf_supVIP,mkdr_supVIP,mka_supVIP,
     &  hka_supVIP,mk2_supVIP,hk2_supVIP,
     &  mkm_supVIP,mkc_supVIP,mkahp_supVIP,
     &  mcat_supVIP,hcat_supVIP,mcal_supVIP,
     &  mar_supVIP)


        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
       do L = 1, num_supVIP     
        distal_axon_supintern   (L + 300     ) = V_supVIP      (59,L)
       end do
  
           call mpi_allgather (distal_axon_supintern,
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_sup_local(1) = 0.d0     
        field_deep_local(1) = 0.d0     
           call mpi_allgather (field_sup_local,     
     &  1              , mpi_double_precision,
     &  field_sup_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_deep_local,     
     &  1              , mpi_double_precision,
     &  field_deep_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
        ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for supVIP 
 
c      ELSE IF (THISNO.EQ.5) THEN
       ELSE IF (nodecell(thisno) .eq. 'LOT         ') THEN
c LOT

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 
          lastcell = num_LOT                             

       IF (mod(O,how_often).eq.0) then
c 1st set LOT synaptic conductances to 0:

          do i = 1, numcomp_LOT
          do j = firstcell, lastcell  
         gAMPA_LOT(i,j)   = 0.d0
         gNMDA_LOT(i,j)   = 0.d0
         gGABA_A_LOT(i,j) = 0.d0
         gGABA_B_LOT(i,j) = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle L2pyr    -> LOT
      do i = 1, num_L2pyr_to_LOT
       j = map_L2pyr_to_LOT(i,L) ! j = presynaptic cell
       k = com_L2pyr_to_LOT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L2pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L2pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L2pyr_to_LOT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_LOT(k,L)  = gAMPA_LOT(k,L) +
     &  gAMPA_L2pyr_to_LOT * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_LOT(k,L) = gNMDA_LOT(k,L) +
     &  gNMDA_L2pyr_to_LOT * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L2pyr_to_LOT
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_LOT(k,L) = gNMDA_LOT(k,L) +
     &  gNMDA_L2pyr_to_LOT * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L2pyr_to_LOT
       if (gNMDA_LOT(k,L).gt.z)
     &  gNMDA_LOT(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle placeholder1     -> LOT
      do i = 1, num_placeholder1_to_LOT
       j = map_placeholder1_to_LOT(i,L) ! j = presynaptic cell
       k = com_placeholder1_to_LOT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder1(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder1(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder1_to_LOT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_LOT(k,L)  = gGABA_A_LOT(k,L) +
     &  gGABA_placeholder1_to_LOT * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle placeholder2     -> LOT
      do i = 1, num_placeholder2_to_LOT
       j = map_placeholder2_to_LOT(i,L) ! j = presynaptic cell
       k = com_placeholder2_to_LOT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder2(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder2(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder2_to_LOT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_LOT(k,L)  = gGABA_A_LOT(k,L) +
     &  gGABA_placeholder2_to_LOT * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle placeholder3      -> LOT
      do i = 1, num_placeholder3_to_LOT
       j = map_placeholder3_to_LOT(i,L) ! j = presynaptic cell
       k = com_placeholder3_to_LOT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder3(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder3(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder3_to_LOT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_LOT(k,L)  = gGABA_A_LOT(k,L) +
     &  gGABA_placeholder3_to_LOT * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle LOT   -> LOT
      do i = 1, num_LOT_to_LOT
       j = map_LOT_to_LOT(i,L) ! j = presynaptic cell
       k = com_LOT_to_LOT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_LOT(j)  ! enumerate presyn. spikes
        presyntime = outtime_LOT(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_LOT_to_LOT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_LOT(k,L)  = gAMPA_LOT(k,L) +
     &  gAMPA_LOT_to_LOT * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_LOT(k,L) = gNMDA_LOT(k,L) +
     &  gNMDA_LOT_to_LOT * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_LOT_to_LOT
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_LOT(k,L) = gNMDA_LOT(k,L) +
     &  gNMDA_LOT_to_LOT * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_LOT_to_LOT
       if (gNMDA_LOT(k,L).gt.z)
     &  gNMDA_LOT(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle semilunar      -> LOT
      do i = 1, num_semilunar_to_LOT
       j = map_semilunar_to_LOT(i,L) ! j = presynaptic cell
       k = com_semilunar_to_LOT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_semilunar(j)  ! enumerate presyn. spikes
        presyntime = outtime_semilunar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_semilunar_to_LOT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_LOT(k,L)  = gAMPA_LOT(k,L) +
     &  gAMPA_semilunar_to_LOT * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_LOT(k,L) = gNMDA_LOT(k,L) +
     &  gNMDA_semilunar_to_LOT * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_semilunar_to_LOT
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_LOT(k,L) = gNMDA_LOT(k,L) +
     &  gNMDA_semilunar_to_LOT * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_semilunar_to_LOT
       if (gNMDA_LOT(k,L).gt.z)
     &  gNMDA_LOT(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle multipolar      -> LOT
      do i = 1, num_multipolar_to_LOT
       j = map_multipolar_to_LOT(i,L) ! j = presynaptic cell
       k = com_multipolar_to_LOT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_multipolar(j)  ! enumerate presyn. spikes
        presyntime = outtime_multipolar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_multipolar_to_LOT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_LOT(k,L)  = gAMPA_LOT(k,L) +
     &  gAMPA_multipolar_to_LOT * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_LOT(k,L) = gNMDA_LOT(k,L) +
     &  gNMDA_multipolar_to_LOT * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_multipolar_to_LOT
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_LOT(k,L) = gNMDA_LOT(k,L) +
     &  gNMDA_multipolar_to_LOT * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_multipolar_to_LOT
       if (gNMDA_LOT(k,L).gt.z)
     &  gNMDA_LOT(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask    -> LOT
      do i = 1, num_deepbask_to_LOT
       j = map_deepbask_to_LOT(i,L) ! j = presynaptic cell
       k = com_deepbask_to_LOT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_LOT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_LOT(k,L)  = gGABA_A_LOT(k,L) +
     &  gGABA_deepbask_to_LOT * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle deepng     -> LOT
      do i = 1, num_deepng_to_LOT
       j = map_deepng_to_LOT(i,L) ! j = presynaptic cell
       k = com_deepng_to_LOT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepng(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_deepng_to_LOT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_LOT(k,L)  = gGABA_A_LOT(k,L) +
     &  gGABA_deepng_to_LOT * z      
! end GABA-A part

        dexparg = delta / tauGABAB_deepng_to_LOT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_B_LOT(k,L)  = gGABA_B_LOT(k,L) +
     &  gGABAB_deepng_to_LOT * z      
! end GABA-A part

c     gGABA_B_LOT(k,L) = gGABA_B_LOT(k,L) +
c    &   gGABAB_deepng_to_LOT * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i



c Handle deepLTS    -> LOT
      do i = 1, num_deepLTS_to_LOT
       j = map_deepLTS_to_LOT(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_LOT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_LOT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_LOT(k,L)  = gGABA_A_LOT(k,L) +
     &  gGABA_deepLTS_to_LOT * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supVIP     -> LOT
      do i = 1, num_supVIP_to_LOT
       j = map_supVIP_to_LOT(i,L) ! j = presynaptic cell
       k = com_supVIP_to_LOT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supVIP(j)  ! enumerate presyn. spikes
        presyntime = outtime_supVIP(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supVIP_to_LOT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_LOT(k,L)  = gGABA_A_LOT(k,L) +
     &  gGABA_supVIP_to_LOT * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle placeholder5         -> LOT
      do i = 1, num_placeholder5_to_LOT
       j = map_placeholder5_to_LOT(i,L) ! j = presynaptic cell
       k = com_placeholder5_to_LOT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder5(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder5(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_placeholder5_to_LOT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_LOT(k,L)  = gAMPA_LOT(k,L) +
     &  gAMPA_placeholder5_to_LOT * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_LOT(k,L) = gNMDA_LOT(k,L) +
     &  gNMDA_placeholder5_to_LOT * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_placeholder5_to_LOT
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_LOT(k,L) = gNMDA_LOT(k,L) +
     &  gNMDA_placeholder5_to_LOT * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_placeholder5_to_LOT 
       if (gNMDA_LOT(k,L).gt.z)
     &  gNMDA_LOT(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle L3pyr   -> LOT
      do i = 1, num_L3pyr_to_LOT
       j = map_L3pyr_to_LOT(i,L) ! j = presynaptic cell
       k = com_L3pyr_to_LOT(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L3pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L3pyr_to_LOT
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_LOT(k,L)  = gAMPA_LOT(k,L) +
     &  gAMPA_L3pyr_to_LOT * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_LOT(k,L) = gNMDA_LOT(k,L) +
     &  gNMDA_L3pyr_to_LOT * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L3pyr_to_LOT
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_LOT(k,L) = gNMDA_LOT(k,L) +
     &  gNMDA_L3pyr_to_LOT * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L3pyr_to_LOT
       if (gNMDA_LOT(k,L).gt.z)
     &  gNMDA_LOT(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of LOT
       ENDIF ! if (mod(O,how_often).eq.0) ...

! Define currents to LOT cells, ectopic spikes,
! tonic synaptic conductances

c       IF (time.ge.400.d0) THEN
c     IF ((time.ge.600.d0).and.(time.le.800.d0)) THEN
      IF ((time.ge.700.d0).and.(time.le.900.d0)) THEN
c     IF  (time.ge.400.d0)                       THEN
c     IF  (time.ge.1000.d0)                       THEN
      if (mod(O,200).eq.0) then
       call durand(seed,num_LOT,ranvec_LOT) 
        do L = firstcell, lastcell  ! Note
         if ((ranvec_LOT(L).gt.0.d0).and.
     &     (ranvec_LOT(L).le.noisepe_LOT)) then
          curr_LOT(59,L) = 0.8d0
c         write(6,2299) 'ectopic', time, L
2299      format (a8,f9.2,i5)
         else
          curr_LOT(59,L) = 0.d0
         endif 
        end do
      endif
        ENDIF


! Call integration routine for LOT cells
       CALL INTEGRATE_LOT (O, time, num_LOT,
     &    V_LOT, curr_LOT,
     &    initialize, firstcell, lastcell,
     & gAMPA_LOT, gNMDA_LOT, gGABA_A_LOT,
     & gGABA_B_LOT, Mg, 
     & gapcon_LOT,totaxgj_LOT,gjtable_LOT, dt,
     &  chi_LOT,mnaf_LOT,mnap_LOT,
     &  hnaf_LOT,mkdr_LOT,mka_LOT,
     &  hka_LOT,mk2_LOT,hk2_LOT,
     &  mkm_LOT,mkc_LOT,mkahp_LOT,
     &  mcat_LOT,hcat_LOT,mcal_LOT,
     &  mar_LOT)


       IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_LOT
       do L = firstcell, lastcell
        distal_axon_LOT (L-firstcell+1) = V_LOT (59,L)
       end do
  
           call mpi_allgather (distal_axon_LOT,
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_sup_local(1) = 0.d0     
        field_deep_local(1) = 0.d0     
           call mpi_allgather (field_sup_local,     
     &  1              , mpi_double_precision,
     &  field_sup_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_deep_local,     
     &  1              , mpi_double_precision,
     &  field_deep_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
           ENDIF !  if (mod(O,how_often).eq.0) ...

! END thisno for LOT

c      else if (1.eq.2) then  ! this should prevent semilunar from
c processing
       ELSE IF (nodecell(thisno) .eq. 'semilunar   ') THEN
c semilunar

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 
          lastcell = num_semilunar                            

         IF (mod(O,how_often).eq.0) then
c 1st set semilunar    synaptic conductances to 0:

          do i = 1, numcomp_semilunar
          do j = firstcell, lastcell
         gAMPA_semilunar(i,j)      = 0.d0
         gNMDA_semilunar(i,j)      = 0.d0
         gGABA_A_semilunar(i,j)    = 0.d0
         gGABA_B_semilunar(i,j)    = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle L2pyr    -> semilunar
      do i = 1, num_L2pyr_to_semilunar   
       j = map_L2pyr_to_semilunar(i,L) ! j = presynaptic cell
       k = com_L2pyr_to_semilunar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L2pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L2pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L2pyr_to_semilunar   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_semilunar(k,L)  = gAMPA_semilunar(k,L) +
     &  gAMPA_L2pyr_to_semilunar * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_semilunar(k,L) = gNMDA_semilunar(k,L) +
     &  gNMDA_L2pyr_to_semilunar * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L2pyr_to_semilunar   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_semilunar(k,L) = gNMDA_semilunar(k,L) +
     &  gNMDA_L2pyr_to_semilunar * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L2pyr_to_semilunar
       if (gNMDA_semilunar(k,L).gt.z)
     &  gNMDA_semilunar(k,L) = z
! end NMDA part

       end do ! m
      end do ! i

c Handle supng      -> semilunar
      do i = 1, num_supng_to_semilunar
       j = map_supng_to_semilunar(i,L) ! j = presynaptic cell
       k = com_supng_to_semilunar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supng(j)  ! enumerate presyn. spikes
        presyntime = outtime_supng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_supng_to_semilunar
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_semilunar(k,L)  = gGABA_A_semilunar(k,L) +
     &  gGABA_supng_to_semilunar * z      
! end GABA-A part

        dexparg = delta / tauGABAB_supng_to_semilunar
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_B_semilunar(k,L)  = gGABA_B_semilunar(k,L) +
     &  gGABAB_supng_to_semilunar * z      
! end GABA-A part

c     gGABA_B_semilunar(k,L) = gGABA_B_semilunar(k,L) +
c    &   gGABAB_supng_to_semilunar * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle placeholder2     -> semilunar
      do i = 1, num_placeholder2_to_semilunar   
       j = map_placeholder2_to_semilunar(i,L) ! j = presynaptic cell
       k = com_placeholder2_to_semilunar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder2(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder2(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder2_to_semilunar   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_semilunar(k,L)  = gGABA_A_semilunar(k,L) +
     &  gGABA_placeholder2_to_semilunar * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle placeholder3      -> semilunar
      do i = 1, num_placeholder3_to_semilunar   
       j = map_placeholder3_to_semilunar(i,L) ! j = presynaptic cell
       k = com_placeholder3_to_semilunar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder3(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder3(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder3_to_semilunar   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_semilunar(k,L)  = gGABA_A_semilunar(k,L) +
     &  gGABA_placeholder3_to_semilunar * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle LOT   -> semilunar
      do i = 1, num_LOT_to_semilunar  
       j = map_LOT_to_semilunar(i,L) ! j = presynaptic cell
       k = com_LOT_to_semilunar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_LOT(j)  ! enumerate presyn. spikes
        presyntime = outtime_LOT(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_LOT_to_semilunar  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_semilunar(k,L)  = gAMPA_semilunar(k,L) +
     &  gAMPA_LOT_to_semilunar * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_semilunar(k,L) = gNMDA_semilunar(k,L) +
     &  gNMDA_LOT_to_semilunar * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_LOT_to_semilunar  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_semilunar(k,L) = gNMDA_semilunar(k,L) +
     &  gNMDA_LOT_to_semilunar * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_LOT_to_semilunar  
       if (gNMDA_semilunar(k,L).gt.z)
     &  gNMDA_semilunar(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle semilunar      -> semilunar
      do i = 1, num_semilunar_to_semilunar  
       j = map_semilunar_to_semilunar(i,L) ! j = presynaptic cell
       k = com_semilunar_to_semilunar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_semilunar(j)  ! enumerate presyn. spikes
        presyntime = outtime_semilunar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_semilunar_to_semilunar  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_semilunar(k,L)  = gAMPA_semilunar(k,L) +
     &  gAMPA_semilunar_to_semilunar * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_semilunar(k,L) = gNMDA_semilunar(k,L) +
     &  gNMDA_semilunar_to_semilunar * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_semilunar_to_semilunar  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_semilunar(k,L) = gNMDA_semilunar(k,L) +
     &  gNMDA_semilunar_to_semilunar * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_semilunar_to_semilunar  
       if (gNMDA_semilunar(k,L).gt.z)
     &  gNMDA_semilunar(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle multipolar      -> semilunar
      do i = 1, num_multipolar_to_semilunar   
       j = map_multipolar_to_semilunar(i,L) ! j = presynaptic cell
       k = com_multipolar_to_semilunar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_multipolar(j)  ! enumerate presyn. spikes
        presyntime = outtime_multipolar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_multipolar_to_semilunar
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_semilunar(k,L)  = gAMPA_semilunar(k,L) +
     &  gAMPA_multipolar_to_semilunar * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_semilunar(k,L) = gNMDA_semilunar(k,L) +
     &  gNMDA_multipolar_to_semilunar * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_multipolar_to_semilunar
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_semilunar(k,L) = gNMDA_semilunar(k,L) +
     &  gNMDA_multipolar_to_semilunar * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_multipolar_to_semilunar   
       if (gNMDA_semilunar(k,L).gt.z)
     &  gNMDA_semilunar(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask    -> semilunar
      do i = 1, num_deepbask_to_semilunar   
       j = map_deepbask_to_semilunar(i,L) ! j = presynaptic cell
       k = com_deepbask_to_semilunar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_semilunar   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_semilunar(k,L)  = gGABA_A_semilunar(k,L) +
     &  gGABA_deepbask_to_semilunar * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepng      -> semilunar
      do i = 1, num_deepng_to_semilunar
       j = map_deepng_to_semilunar(i,L) ! j = presynaptic cell
       k = com_deepng_to_semilunar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepng(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_deepng_to_semilunar
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_semilunar(k,L)  = gGABA_A_semilunar(k,L) +
     &  gGABA_deepng_to_semilunar * z      
! end GABA-A part

        dexparg = delta / tauGABAB_deepng_to_semilunar
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_B_semilunar(k,L)  = gGABA_B_semilunar(k,L) +
     &  gGABAB_deepng_to_semilunar * z      

c     gGABA_B_semilunar(k,L) = gGABA_B_semilunar(k,L) +
c    &   gGABAB_deepng_to_semilunar * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle deepLTS    -> semilunar
      do i = 1, num_deepLTS_to_semilunar   
       j = map_deepLTS_to_semilunar(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_semilunar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_semilunar   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_semilunar(k,L)  = gGABA_A_semilunar(k,L) +
     &  gGABA_deepLTS_to_semilunar * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supVIP     -> semilunar
      do i = 1, num_supVIP_to_semilunar   
       j = map_supVIP_to_semilunar(i,L) ! j = presynaptic cell
       k = com_supVIP_to_semilunar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supVIP(j)  ! enumerate presyn. spikes
        presyntime = outtime_supVIP(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta)

! GABA-A part
        dexparg = delta / tauGABA_supVIP_to_semilunar   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_semilunar(k,L)  = gGABA_A_semilunar(k,L) +
     &  gGABA_supVIP_to_semilunar * z      
! end GABA-A part

c  k0 must be properly defined
      gGABA_B_semilunar  (k,L) = gGABA_B_semilunar  (k,L) +
     &   gGABAB_supVIP_to_semilunar   * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle placeholder5         -> semilunar
      do i = 1, num_placeholder5_to_semilunar
       j = map_placeholder5_to_semilunar(i,L) ! j = presynaptic cell
       k = com_placeholder5_to_semilunar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder5(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder5(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_placeholder5_to_semilunar
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_semilunar(k,L)  = gAMPA_semilunar(k,L) +
     &  gAMPA_placeholder5_to_semilunar * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_semilunar(k,L) = gNMDA_semilunar(k,L) +
     &  gNMDA_placeholder5_to_semilunar * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_placeholder5_to_semilunar
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_semilunar(k,L) = gNMDA_semilunar(k,L) +
     &  gNMDA_placeholder5_to_semilunar * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_placeholder5_to_semilunar 
       if (gNMDA_semilunar(k,L).gt.z)
     &  gNMDA_semilunar(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle L3pyr   -> semilunar
      do i = 1, num_L3pyr_to_semilunar
       j = map_L3pyr_to_semilunar(i,L) ! j = presynaptic cell
       k = com_L3pyr_to_semilunar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L3pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L3pyr_to_semilunar   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_semilunar(k,L)  = gAMPA_semilunar(k,L) +
     &  gAMPA_L3pyr_to_semilunar * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_semilunar(k,L) = gNMDA_semilunar(k,L) +
     &  gNMDA_L3pyr_to_semilunar * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L3pyr_to_semilunar   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_semilunar(k,L) = gNMDA_semilunar(k,L) +
     &  gNMDA_L3pyr_to_semilunar * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L3pyr_to_semilunar   
       if (gNMDA_semilunar(k,L).gt.z)
     &  gNMDA_semilunar(k,L) = z
! end NMDA part

       end do ! m
      end do ! i

         end do
c End enumeration of semilunar   
         ENDIF  ! if (mod(O,how_often).eq.0) ....

! Define currents to semilunar    cells, ectopic spikes,
! tonic synaptic conductances

      if (mod(O,200).eq.0) then
       call durand(seed,num_semilunar  ,ranvec_semilunar  ) 
        do L = firstcell, lastcell
         if ((ranvec_semilunar  (L).gt.0.d0).and.
     &     (ranvec_semilunar  (L).le.noisepe_semilunar  )) then
          curr_semilunar  (60,L) = 0.4d0
         else
          curr_semilunar  (60,L) = 0.d0
         endif 
        end do
      endif

! Call integration routine for semilunar    cells
       CALL INTEGRATE_semilunar (O, time, num_semilunar,
     &    V_semilunar, curr_semilunar,
     &    initialize, firstcell, lastcell,
     & gAMPA_semilunar, gNMDA_semilunar, gGABA_A_semilunar,
     & gGABA_B_semilunar, Mg, 
     & gapcon_semilunar  ,totaxgj_semilunar   ,gjtable_semilunar, dt,
     &  chi_semilunar,mnaf_semilunar,mnap_semilunar,
     &  hnaf_semilunar,mkdr_semilunar,mka_semilunar,
     &  hka_semilunar,mk2_semilunar,hk2_semilunar,
     &  mkm_semilunar,mkc_semilunar,mkahp_semilunar,
     &  mcat_semilunar,hcat_semilunar,mcal_semilunar,
     &  mar_semilunar,field_sup ,field_deep, rel_axonshift_semilunar)
  

        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_semilunar   
       do L = firstcell, lastcell
        distal_axon_semilunar    (L-firstcell+1) = V_semilunar    (72,L)
       end do
  
           call mpi_allgather (distal_axon_semilunar,  
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_sup_local(1) = field_sup
        field_deep_local(1) = field_deep

           call mpi_allgather (field_sup_local,     
     &  1              , mpi_double_precision,
     &  field_sup_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

           call mpi_allgather (field_deep_local,     
     &  1              , mpi_double_precision,
     &  field_deep_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
4000         continue
  
           ENDIF  ! if (mod(O,how_often).eq.0) ...

! END thisno for semilunar

c      ELSE IF (THISNO.EQ.7) THEN
       ELSE IF (nodecell(thisno) .eq. 'multipolar') THEN
c multipolar

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 
          lastcell = num_multipolar                             

         IF (mod(O,how_often).eq.0) then
c 1st set multipolar    synaptic conductances to 0:

          do i = 1, numcomp_multipolar
          do j = firstcell, lastcell
         gAMPA_multipolar(i,j)      = 0.d0 
         gNMDA_multipolar(i,j)      = 0.d0
         gGABA_A_multipolar(i,j)    = 0.d0
c        gGABA_B_multipolar(i,j)    = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle L2pyr    -> multipolar
      do i = 1, num_L2pyr_to_multipolar   
       j = map_L2pyr_to_multipolar(i,L) ! j = presynaptic cell
       k = com_L2pyr_to_multipolar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L2pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L2pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L2pyr_to_multipolar   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_multipolar(k,L)  = gAMPA_multipolar(k,L) +
     &  gAMPA_L2pyr_to_multipolar * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_multipolar(k,L) = gNMDA_multipolar(k,L) +
     &  gNMDA_L2pyr_to_multipolar * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L2pyr_to_multipolar   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_multipolar(k,L) = gNMDA_multipolar(k,L) +
     &  gNMDA_L2pyr_to_multipolar * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L2pyr_to_multipolar
       if (gNMDA_multipolar(k,L).gt.z)
     &  gNMDA_multipolar(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supng      -> multipolar
      do i = 1, num_supng_to_multipolar
       j = map_supng_to_multipolar(i,L) ! j = presynaptic cell
       k = com_supng_to_multipolar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supng(j)  ! enumerate presyn. spikes
        presyntime = outtime_supng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_supng_to_multipolar
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_multipolar(k,L)  = gGABA_A_multipolar(k,L) +
     &  gGABA_supng_to_multipolar * z      
! end GABA-A part

        dexparg = delta / tauGABAB_supng_to_multipolar
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

c     gGABA_B_multipolar(k,L)  = gGABA_B_multipolar(k,L) +
c    &  gGABAB_supng_to_multipolar * z      

c     gGABA_B_multipolar(k,L) = gGABA_B_multipolar(k,L) +
c    &   gGABAB_supng_to_multipolar * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle placeholder2     -> multipolar
      do i = 1, num_placeholder2_to_multipolar   
       j = map_placeholder2_to_multipolar(i,L) ! j = presynaptic cell
       k = com_placeholder2_to_multipolar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder2(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder2(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder2_to_multipolar   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_multipolar(k,L)  = gGABA_A_multipolar(k,L) +
     &  gGABA_placeholder2_to_multipolar * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle placeholder3      -> multipolar
      do i = 1, num_placeholder3_to_multipolar   
       j = map_placeholder3_to_multipolar(i,L) ! j = presynaptic cell
       k = com_placeholder3_to_multipolar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder3(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder3(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder3_to_multipolar   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_multipolar(k,L)  = gGABA_A_multipolar(k,L) +
     &  gGABA_placeholder3_to_multipolar * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle LOT   -> multipolar
      do i = 1, num_LOT_to_multipolar  
       j = map_LOT_to_multipolar(i,L) ! j = presynaptic cell
       k = com_LOT_to_multipolar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_LOT(j)  ! enumerate presyn. spikes
        presyntime = outtime_LOT(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_LOT_to_multipolar  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_multipolar(k,L)  = gAMPA_multipolar(k,L) +
     &  gAMPA_LOT_to_multipolar * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_multipolar(k,L) = gNMDA_multipolar(k,L) +
     &  gNMDA_LOT_to_multipolar * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_LOT_to_multipolar  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_multipolar(k,L) = gNMDA_multipolar(k,L) +
     &  gNMDA_LOT_to_multipolar * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_LOT_to_multipolar  
       if (gNMDA_multipolar(k,L).gt.z)
     &  gNMDA_multipolar(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle semilunar      -> multipolar
      do i = 1, num_semilunar_to_multipolar  
       j = map_semilunar_to_multipolar(i,L) ! j = presynaptic cell
       k = com_semilunar_to_multipolar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_semilunar(j)  ! enumerate presyn. spikes
        presyntime = outtime_semilunar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_semilunar_to_multipolar  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_multipolar(k,L)  = gAMPA_multipolar(k,L) +
     &  gAMPA_semilunar_to_multipolar * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_multipolar(k,L) = gNMDA_multipolar(k,L) +
     &  gNMDA_semilunar_to_multipolar * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_semilunar_to_multipolar  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_multipolar(k,L) = gNMDA_multipolar(k,L) +
     &  gNMDA_semilunar_to_multipolar * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_semilunar_to_multipolar  
       if (gNMDA_multipolar(k,L).gt.z)
     &  gNMDA_multipolar(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle multipolar      -> multipolar
      do i = 1, num_multipolar_to_multipolar  
       j = map_multipolar_to_multipolar(i,L) ! j = presynaptic cell
       k = com_multipolar_to_multipolar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_multipolar(j)  ! enumerate presyn. spikes
        presyntime = outtime_multipolar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_multipolar_to_multipolar  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_multipolar(k,L)  = gAMPA_multipolar(k,L) +
     &  gAMPA_multipolar_to_multipolar * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_multipolar(k,L) = gNMDA_multipolar(k,L) +
     &  gNMDA_multipolar_to_multipolar * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_multipolar_to_multipolar  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_multipolar(k,L) = gNMDA_multipolar(k,L) +
     &  gNMDA_multipolar_to_multipolar * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_multipolar_to_multipolar  
       if (gNMDA_multipolar(k,L).gt.z)
     &  gNMDA_multipolar(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask    -> multipolar
      do i = 1, num_deepbask_to_multipolar   
       j = map_deepbask_to_multipolar(i,L) ! j = presynaptic cell
       k = com_deepbask_to_multipolar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_multipolar   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_multipolar(k,L)  = gGABA_A_multipolar(k,L) +
     &  gGABA_deepbask_to_multipolar * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepng      -> multipolar
      do i = 1, num_deepng_to_multipolar
       j = map_deepng_to_multipolar(i,L) ! j = presynaptic cell
       k = com_deepng_to_multipolar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepng(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_deepng_to_multipolar
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_multipolar(k,L)  = gGABA_A_multipolar(k,L) +
     &  gGABA_deepng_to_multipolar * z      
! end GABA-A part

c       dexparg = delta / tauGABAB_deepng_to_multipolar
c note that dexparg = MINUS the actual arg. to dexp
c        if (dexparg.le.5.d0) then
c         z = dexptablesmall (int(dexparg*1000.d0))
c        else if (dexparg.le.100.d0) then
c         z = dexptablebig (int(dexparg*10.d0))
c        else
c         z = 0.d0
c        endif

c     gGABA_B_multipolar(k,L)  = gGABA_B_multipolar(k,L) +
c    &  gGABAB_deepng_to_multipolar * z      

c     gGABA_B_multipolar(k,L) = gGABA_B_multipolar(k,L) +
c    &   gGABAB_deepng_to_multipolar * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle deepLTS    -> multipolar
      do i = 1, num_deepLTS_to_multipolar   
       j = map_deepLTS_to_multipolar(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_multipolar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_multipolar   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_multipolar(k,L)  = gGABA_A_multipolar(k,L) +
     &  gGABA_deepLTS_to_multipolar * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supVIP     -> multipolar
      do i = 1, num_supVIP_to_multipolar   
       j = map_supVIP_to_multipolar(i,L) ! j = presynaptic cell
       k = com_supVIP_to_multipolar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supVIP(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta)

! GABA-A part
        dexparg = delta / tauGABA_supVIP_to_multipolar   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_multipolar(k,L)  = gGABA_A_multipolar(k,L) +
     &  gGABA_supVIP_to_multipolar * z      
! end GABA-A part

c  k0 must be properly defined
c     gGABA_B_multipolar(k,L) = gGABA_B_multipolar(k,L) +
c    &   gGABAB_supVIP_to_multipolar * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle placeholder5         -> multipolar
      do i = 1, num_placeholder5_to_multipolar
       j = map_placeholder5_to_multipolar(i,L) ! j = presynaptic cell
       k = com_placeholder5_to_multipolar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder5(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder5(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_placeholder5_to_multipolar
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_multipolar(k,L)  = gAMPA_multipolar(k,L) +
     &  gAMPA_placeholder5_to_multipolar * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_multipolar(k,L) = gNMDA_multipolar(k,L) +
     &  gNMDA_placeholder5_to_multipolar * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_placeholder5_to_multipolar
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_multipolar(k,L) = gNMDA_multipolar(k,L) +
     &  gNMDA_placeholder5_to_multipolar * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_placeholder5_to_multipolar 
       if (gNMDA_multipolar(k,L).gt.z)
     &  gNMDA_multipolar(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle L3pyr   -> multipolar
      do i = 1, num_L3pyr_to_multipolar  
       j = map_L3pyr_to_multipolar(i,L) ! j = presynaptic cell
       k = com_L3pyr_to_multipolar(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L3pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L3pyr_to_multipolar  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_multipolar(k,L)  = gAMPA_multipolar(k,L) +
     &  gAMPA_L3pyr_to_multipolar * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_multipolar(k,L) = gNMDA_multipolar(k,L) +
     &  gNMDA_L3pyr_to_multipolar * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L3pyr_to_multipolar  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_multipolar(k,L) = gNMDA_multipolar(k,L) +
     &  gNMDA_L3pyr_to_multipolar * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L3pyr_to_multipolar  
       if (gNMDA_multipolar(k,L).gt.z)
     &  gNMDA_multipolar(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of multipolar   
        ENDIF  ! if (mod(O,how_often).eq.0) ...

! Define currents to multipolar    cells, ectopic spikes,
! tonic synaptic conductances

      if (mod(O,200).eq.0) then
       call durand(seed,num_multipolar  ,ranvec_multipolar  ) 
        do L = firstcell, lastcell
         if ((ranvec_multipolar  (L).gt.0.d0).and.
     &     (ranvec_multipolar  (L).le.noisepe_multipolar  )) then
          curr_multipolar  (58,L) = 0.4d0
         else
          curr_multipolar  (58,L) = 0.d0
         endif 
        end do
      endif

! Call integration routine for multipolar    cells
       CALL INTEGRATE_multipolar (O, time, num_multipolar,
     &    V_multipolar, curr_multipolar,
     & initialize, firstcell, lastcell,
     & gAMPA_multipolar, gNMDA_multipolar, gGABA_A_multipolar,
c    & gGABA_B_multipolar, Mg, 
     &                     Mg, 
     & gapcon_multipolar,totaxgj_multipolar,gjtable_multipolar,dt,
     &  chi_multipolar,mnaf_multipolar,mnap_multipolar,
     &  hnaf_multipolar,mkdr_multipolar,mka_multipolar,
     &  hka_multipolar,mk2_multipolar,hk2_multipolar,
     &  mkm_multipolar,mkc_multipolar,mkahp_multipolar,
     &  mcat_multipolar,hcat_multipolar,mcal_multipolar,
     &  mar_multipolar )
c    &  mar_multipolar,field_sup       ,field_deep       )

  

       IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_multipolar   
       do L = firstcell, lastcell
        distal_axon_multipolar (L-firstcell+1) = V_multipolar(59,L)
       end do
  
           call mpi_allgather (distal_axon_multipolar,  
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

c       field_sup_local(1) = field_sup
        field_sup_local(1) = 0.d0      
c       field_deep_local(1) = field_deep
        field_deep_local(1) = 0.d0          
           call mpi_allgather (field_sup_local,     
     &  1              , mpi_double_precision,
     &  field_sup_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_deep_local,     
     &  1              , mpi_double_precision,
     &  field_deep_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
         ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for multipolar

c      ELSE IF (THISNO.EQ.8) THEN
       ELSE IF (nodecell(thisno) .eq. 'L3pyr       ') THEN
c L3pyr

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 
          lastcell = num_L3pyr                            

         IF (mod(O,how_often).eq.0) then
c 1st set L3pyr synaptic conductances to 0:

          do i = 1, numcomp_L3pyr
          do j = firstcell, lastcell
         gAMPA_L3pyr(i,j)   = 0.d0 
         gNMDA_L3pyr(i,j)   = 0.d0 
         gGABA_A_L3pyr(i,j) = 0.d0
         gGABA_B_L3pyr(i,j) = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle L2pyr   -> L3pyr
      do i = 1, num_L2pyr_to_L3pyr   
       j = map_L2pyr_to_L3pyr(i,L) ! j = presynaptic cell
       k = com_L2pyr_to_L3pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L2pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L2pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L2pyr_to_L3pyr   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_L3pyr(k,L)  = gAMPA_L3pyr(k,L) +
     &  gAMPA_L2pyr_to_L3pyr * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_L3pyr(k,L) = gNMDA_L3pyr(k,L) +
     &  gNMDA_L2pyr_to_L3pyr * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L2pyr_to_L3pyr   
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_L3pyr(k,L) = gNMDA_L3pyr(k,L) +
     &  gNMDA_L2pyr_to_L3pyr * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L2pyr_to_L3pyr
       if (gNMDA_L3pyr(k,L).gt.z)
     &  gNMDA_L3pyr(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle supng      -> L3pyr
      do i = 1, num_supng_to_L3pyr
       j = map_supng_to_L3pyr(i,L) ! j = presynaptic cell
       k = com_supng_to_L3pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supng(j)  ! enumerate presyn. spikes
        presyntime = outtime_supng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_supng_to_L3pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L3pyr(k,L)  = gGABA_A_L3pyr(k,L) +
     &  gGABA_supng_to_L3pyr * z      
! end GABA-A part

        dexparg = delta / tauGABAB_supng_to_L3pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_B_L3pyr(k,L)  = gGABA_B_L3pyr(k,L) +
     &  gGABAB_supng_to_L3pyr * z      

c     gGABA_B_L3pyr(k,L) = gGABA_B_L3pyr(k,L) +
c    &   gGABAB_supng_to_L3pyr * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle placeholder2    -> L3pyr
      do i = 1, num_placeholder2_to_L3pyr   
       j = map_placeholder2_to_L3pyr(i,L) ! j = presynaptic cell
       k = com_placeholder2_to_L3pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder2(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder2(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder2_to_L3pyr   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L3pyr(k,L)  = gGABA_A_L3pyr(k,L) +
     &  gGABA_placeholder2_to_L3pyr * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle placeholder3     -> L3pyr
      do i = 1, num_placeholder3_to_L3pyr   
       j = map_placeholder3_to_L3pyr(i,L) ! j = presynaptic cell
       k = com_placeholder3_to_L3pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder3(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder3(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder3_to_L3pyr   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L3pyr(k,L)  = gGABA_A_L3pyr(k,L) +
     &  gGABA_placeholder3_to_L3pyr * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle LOT  -> L3pyr
      do i = 1, num_LOT_to_L3pyr  
       j = map_LOT_to_L3pyr(i,L) ! j = presynaptic cell
       k = com_LOT_to_L3pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_LOT(j)  ! enumerate presyn. spikes
        presyntime = outtime_LOT(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_LOT_to_L3pyr  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_L3pyr(k,L)  = gAMPA_L3pyr(k,L) +
     &  gAMPA_LOT_to_L3pyr * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_L3pyr(k,L) = gNMDA_L3pyr(k,L) +
     &  gNMDA_LOT_to_L3pyr * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_LOT_to_L3pyr  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_L3pyr(k,L) = gNMDA_L3pyr(k,L) +
     &  gNMDA_LOT_to_L3pyr * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_LOT_to_L3pyr  
       if (gNMDA_L3pyr(k,L).gt.z)
     &  gNMDA_L3pyr(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle semilunar     -> L3pyr
      do i = 1, num_semilunar_to_L3pyr  
       j = map_semilunar_to_L3pyr(i,L) ! j = presynaptic cell
       k = com_semilunar_to_L3pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_semilunar(j)  ! enumerate presyn. spikes
        presyntime = outtime_semilunar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_semilunar_to_L3pyr  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_L3pyr(k,L)  = gAMPA_L3pyr(k,L) +
     &  gAMPA_semilunar_to_L3pyr * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_L3pyr(k,L) = gNMDA_L3pyr(k,L) +
     &  gNMDA_semilunar_to_L3pyr * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_semilunar_to_L3pyr  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_L3pyr(k,L) = gNMDA_L3pyr(k,L) +
     &  gNMDA_semilunar_to_L3pyr * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_semilunar_to_L3pyr  
       if (gNMDA_L3pyr(k,L).gt.z)
     &  gNMDA_L3pyr(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle multipolar     -> L3pyr
      do i = 1, num_multipolar_to_L3pyr  
       j = map_multipolar_to_L3pyr(i,L) ! j = presynaptic cell
       k = com_multipolar_to_L3pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_multipolar(j)  ! enumerate presyn. spikes
        presyntime = outtime_multipolar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_multipolar_to_L3pyr  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_L3pyr(k,L)  = gAMPA_L3pyr(k,L) +
     &  gAMPA_multipolar_to_L3pyr * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_L3pyr(k,L) = gNMDA_L3pyr(k,L) +
     &  gNMDA_multipolar_to_L3pyr * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_multipolar_to_L3pyr  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_L3pyr(k,L) = gNMDA_L3pyr(k,L) +
     &  gNMDA_multipolar_to_L3pyr * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_multipolar_to_L3pyr  
       if (gNMDA_L3pyr(k,L).gt.z)
     &  gNMDA_L3pyr(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask   -> L3pyr
      do i = 1, num_deepbask_to_L3pyr   
       j = map_deepbask_to_L3pyr(i,L) ! j = presynaptic cell
       k = com_deepbask_to_L3pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_L3pyr   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L3pyr(k,L)  = gGABA_A_L3pyr(k,L) +
     &  gGABA_deepbask_to_L3pyr * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle deepng      -> L3pyr
      do i = 1, num_deepng_to_L3pyr
       j = map_deepng_to_L3pyr(i,L) ! j = presynaptic cell
       k = com_deepng_to_L3pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepng(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepng(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part AND GABA-B part
        dexparg = delta / tauGABA_deepng_to_L3pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L3pyr(k,L)  = gGABA_A_L3pyr(k,L) +
     &  gGABA_deepng_to_L3pyr * z      
! end GABA-A part

        dexparg = delta / tauGABAB_deepng_to_L3pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_B_L3pyr(k,L)  = gGABA_B_L3pyr(k,L) +
     &  gGABAB_deepng_to_L3pyr * z      

c     gGABA_B_L3pyr(k,L) = gGABA_B_L3pyr(k,L) +
c    &   gGABAB_deepng_to_L3pyr * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle deepLTS   -> L3pyr
      do i = 1, num_deepLTS_to_L3pyr   
       j = map_deepLTS_to_L3pyr(i,L) ! j = presynaptic cell
       k = com_deepLTS_to_L3pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepLTS(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepLTS(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepLTS_to_L3pyr   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L3pyr(k,L)  = gGABA_A_L3pyr(k,L) +
     &  gGABA_deepLTS_to_L3pyr * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supVIP    -> L3pyr
      do i = 1, num_supVIP_to_L3pyr   
       j = map_supVIP_to_L3pyr(i,L) ! j = presynaptic cell
       k = com_supVIP_to_L3pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supVIP(j)  ! enumerate presyn. spikes
        presyntime = outtime_supVIP(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta)

! GABA-A part
        dexparg = delta / tauGABA_supVIP_to_L3pyr   
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_L3pyr(k,L)  = gGABA_A_L3pyr(k,L) +
     &  gGABA_supVIP_to_L3pyr * z      
! end GABA-A part

c  k0 must be properly defined
      gGABA_B_L3pyr(k,L) = gGABA_B_L3pyr(k,L) +
     &   gGABAB_supVIP_to_L3pyr * otis_table(k0)
! end GABA-B part

       end do ! m
      end do ! i


c Handle placeholder5        -> L3pyr
      do i = 1, num_placeholder5_to_L3pyr
       j = map_placeholder5_to_L3pyr(i,L) ! j = presynaptic cell
       k = com_placeholder5_to_L3pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder5(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder5(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_placeholder5_to_L3pyr
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_L3pyr(k,L)  = gAMPA_L3pyr(k,L) +
     &  gAMPA_placeholder5_to_L3pyr * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_L3pyr(k,L) = gNMDA_L3pyr(k,L) +
     &  gNMDA_placeholder5_to_L3pyr * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_placeholder5_to_L3pyr
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_L3pyr(k,L) = gNMDA_L3pyr(k,L) +
     &  gNMDA_placeholder5_to_L3pyr * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_placeholder5_to_L3pyr 
       if (gNMDA_L3pyr(k,L).gt.z)
     &  gNMDA_L3pyr(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle L3pyr  -> L3pyr
      do i = 1, num_L3pyr_to_L3pyr  
       j = map_L3pyr_to_L3pyr(i,L) ! j = presynaptic cell
       k = com_L3pyr_to_L3pyr(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L3pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L3pyr_to_L3pyr  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_L3pyr(k,L)  = gAMPA_L3pyr(k,L) +
     &  gAMPA_L3pyr_to_L3pyr * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_L3pyr(k,L) = gNMDA_L3pyr(k,L) +
     &  gNMDA_L3pyr_to_L3pyr * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L3pyr_to_L3pyr  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_L3pyr(k,L) = gNMDA_L3pyr(k,L) +
     &  gNMDA_L3pyr_to_L3pyr * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L3pyr_to_L3pyr  
       if (gNMDA_L3pyr(k,L).gt.z)
     &  gNMDA_L3pyr(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of L3pyr   
          ENDIF  ! if (mod(O,how_often).eq.0) ...

! Define currents to L3pyr    cells, ectopic spikes,
! tonic synaptic conductances

      if (mod(O,200).eq.0) then
       call durand(seed,num_L3pyr  ,ranvec_L3pyr  ) 
        do L = firstcell, lastcell
         if ((ranvec_L3pyr  (L).gt.0.d0).and.
     &     (ranvec_L3pyr  (L).le.noisepe_L3pyr  )) then
          curr_L3pyr  (48,L) = 0.4d0
         else
          curr_L3pyr  (48,L) = 0.d0
         endif 
        end do
      endif

! Call integration routine for L3pyr    cells
       CALL INTEGRATE_L3pyr (O, time, num_L3pyr,
     &    V_L3pyr, curr_L3pyr,
     &    initialize, firstcell, lastcell,
     & gAMPA_L3pyr, gNMDA_L3pyr, gGABA_A_L3pyr,
     & gGABA_B_L3pyr, Mg, 
     & gapcon_L3pyr  ,totaxgj_L3pyr   ,gjtable_L3pyr, dt,
     &  chi_L3pyr,mnaf_L3pyr,mnap_L3pyr,
     &  hnaf_L3pyr,mkdr_L3pyr,mka_L3pyr,
     &  hka_L3pyr,mk2_L3pyr,hk2_L3pyr,
     &  mkm_L3pyr,mkc_L3pyr,mkahp_L3pyr,
     &  mcat_L3pyr,hcat_L3pyr,mcal_L3pyr,
     &  mar_L3pyr,field_sup ,field_deep, rel_axonshift_L3pyr)


        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_L3pyr   
       do L = firstcell, lastcell
        distal_axon_L3pyr    (L-firstcell+1) = V_L3pyr    (72,L)
       end do
  
           call mpi_allgather (distal_axon_L3pyr,
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_sup_local(1) = field_sup
        field_deep_local(1) = field_deep
           call mpi_allgather (field_sup_local,     
     &  1              , mpi_double_precision,
     &  field_sup_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_deep_local,     
     &  1              , mpi_double_precision,
     &  field_deep_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
         ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for L3pyr

c      ELSE IF (THISNO.EQ.9) THEN
c      ELSE IF (nodecell(thisno) .eq. 'deepbask ') THEN
       ELSE IF (nodecell(thisno) .eq. 'deepintern  ') THEN
c deepbask

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 
          lastcell = num_deepbask                              

          IF (mod(O,how_often).eq.0) then
c 1st set deepbask  synaptic conductances to 0:

          do i = 1, numcomp_deepbask
          do j = firstcell, lastcell
         gAMPA_deepbask(i,j)    = 0.d0
         gNMDA_deepbask(i,j)    = 0.d0
         gGABA_A_deepbask(i,j)  = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle L2pyr   -> deepbask
      do i = 1, num_L2pyr_to_deepbask  
       j = map_L2pyr_to_deepbask(i,L) ! j = presynaptic cell
       k = com_L2pyr_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L2pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L2pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L2pyr_to_deepbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepbask(k,L)  = gAMPA_deepbask(k,L) +
     &  gAMPA_L2pyr_to_deepbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_L2pyr_to_deepbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L2pyr_to_deepbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_L2pyr_to_deepbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L2pyr_to_deepbask  
       if (gNMDA_deepbask(k,L).gt.z)
     &  gNMDA_deepbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle placeholder3     -> deepbask
      do i = 1, num_placeholder3_to_deepbask    
       j = map_placeholder3_to_deepbask(i,L) ! j = presynaptic cell
       k = com_placeholder3_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder3(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder3(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder3_to_deepbask    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepbask(k,L)  = gGABA_A_deepbask(k,L) +
     &  gGABA_placeholder3_to_deepbask * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle LOT  -> deepbask
      do i = 1, num_LOT_to_deepbask   
       j = map_LOT_to_deepbask(i,L) ! j = presynaptic cell
       k = com_LOT_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_LOT(j)  ! enumerate presyn. spikes
        presyntime = outtime_LOT(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_LOT_to_deepbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepbask(k,L)  = gAMPA_deepbask(k,L) +
     &  gAMPA_LOT_to_deepbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_LOT_to_deepbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_LOT_to_deepbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_LOT_to_deepbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_LOT_to_deepbask  
       if (gNMDA_deepbask(k,L).gt.z)
     &  gNMDA_deepbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle semilunar     -> deepbask
      do i = 1, num_semilunar_to_deepbask   
       j = map_semilunar_to_deepbask(i,L) ! j = presynaptic cell
       k = com_semilunar_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_semilunar(j)  ! enumerate presyn. spikes
        presyntime = outtime_semilunar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_semilunar_to_deepbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepbask(k,L)  = gAMPA_deepbask(k,L) +
     &  gAMPA_semilunar_to_deepbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_semilunar_to_deepbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_semilunar_to_deepbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_semilunar_to_deepbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_semilunar_to_deepbask  
       if (gNMDA_deepbask(k,L).gt.z)
     &  gNMDA_deepbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle multipolar     -> deepbask
      do i = 1, num_multipolar_to_deepbask   
       j = map_multipolar_to_deepbask(i,L) ! j = presynaptic cell
       k = com_multipolar_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_multipolar(j)  ! enumerate presyn. spikes
        presyntime = outtime_multipolar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_multipolar_to_deepbask  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepbask(k,L)  = gAMPA_deepbask(k,L) +
     &  gAMPA_multipolar_to_deepbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_multipolar_to_deepbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_multipolar_to_deepbask  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_multipolar_to_deepbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_multipolar_to_deepbask  
       if (gNMDA_deepbask(k,L).gt.z)
     &  gNMDA_deepbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask   -> deepbask
      do i = 1, num_deepbask_to_deepbask    
       j = map_deepbask_to_deepbask(i,L) ! j = presynaptic cell
       k = com_deepbask_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_deepbask    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepbask(k,L)  = gGABA_A_deepbask(k,L) +
     &  gGABA_deepbask_to_deepbask * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle deepng     -> deepbask
      do i = 1, num_deepng_to_deepbask    
       j = map_deepng_to_deepbask(i,L) ! j = presynaptic cell
       k = com_deepng_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepng(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepng(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepng_to_deepbask    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepbask(k,L)  = gGABA_A_deepbask(k,L) +
     &  gGABA_deepng_to_deepbask * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supVIP    -> deepbask
      do i = 1, num_supVIP_to_deepbask    
       j = map_supVIP_to_deepbask(i,L) ! j = presynaptic cell
       k = com_supVIP_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supVIP(j)  ! enumerate presyn. spikes
        presyntime = outtime_supVIP(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supVIP_to_deepbask    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepbask(k,L)  = gGABA_A_deepbask(k,L) +
     &  gGABA_supVIP_to_deepbask * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle placeholder5        -> deepbask
      do i = 1, num_placeholder5_to_deepbask 
       j = map_placeholder5_to_deepbask(i,L) ! j = presynaptic cell
       k = com_placeholder5_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder5(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder5(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_placeholder5_to_deepbask 
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepbask(k,L)  = gAMPA_deepbask(k,L) +
     &  gAMPA_placeholder5_to_deepbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_placeholder5_to_deepbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_placeholder5_to_deepbask 
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_placeholder5_to_deepbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_placeholder5_to_deepbask  
       if (gNMDA_deepbask(k,L).gt.z)
     &  gNMDA_deepbask(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle L3pyr  -> deepbask
      do i = 1, num_L3pyr_to_deepbask
       j = map_L3pyr_to_deepbask(i,L) ! j = presynaptic cell
       k = com_L3pyr_to_deepbask(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L3pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L3pyr_to_deepbask
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepbask(k,L)  = gAMPA_deepbask(k,L) +
     &  gAMPA_L3pyr_to_deepbask * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_L3pyr_to_deepbask * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L3pyr_to_deepbask
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepbask(k,L) = gNMDA_deepbask(k,L) +
     &  gNMDA_L3pyr_to_deepbask * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L3pyr_to_deepbask
       if (gNMDA_deepbask(k,L).gt.z)
     &  gNMDA_deepbask(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of deepbask    
         ENDIF ! if (mod(O,how_often).eq.0) ...

! Define currents to deepbask     cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for deepbask     cells
       CALL INTEGRATE_deepbask  (O, time, num_deepbask ,
     &    V_deepbask , curr_deepbask ,
     & initialize, firstcell, lastcell,
     & gAMPA_deepbask, gNMDA_deepbask, gGABA_A_deepbask,
     & Mg, 
     & gapcon_deepbask  ,totSDgj_deepbask   ,gjtable_deepbask, dt,
     &  chi_deepbask,mnaf_deepbask,mnap_deepbask,
     &  hnaf_deepbask,mkdr_deepbask,mka_deepbask,
     &  hka_deepbask,mk2_deepbask,hk2_deepbask,
     &  mkm_deepbask,mkc_deepbask,mkahp_deepbask,
     &  mcat_deepbask,hcat_deepbask,mcal_deepbask,
     &  mar_deepbask)


        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
       do L = 1, num_deepbask    
c      do L = firstcell, lastcell
        distal_axon_deepintern   (L            ) = V_deepbask     (59,L)
       end do
  
c          call mpi_allgather (distal_axon_deepbask,
c    &  maxcellspernode, mpi_double_precision,
c    &  distal_axon_global,maxcellspernode,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)

        field_sup_local(1) = 0.d0     
        field_deep_local(1) = 0.d0     
c          call mpi_allgather (field_sup_local,     
c    &  1              , mpi_double_precision,
c    &  field_sup_global  , 1             ,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)
c          call mpi_allgather (field_deep_local,     
c    &  1              , mpi_double_precision,
c    &  field_deep_global  , 1             ,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)
  
           ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for deepbask

c      ELSE IF (nodecell(thisno) .eq. 'deepng   ') THEN
c deepng  

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 
          lastcell = num_deepng                            

          IF (mod(O,how_often).eq.0) then
c 1st set deepng    synaptic conductances to 0:

          do i = 1, numcomp_deepng  
          do j = firstcell, lastcell
         gAMPA_deepng  (i,j)    = 0.d0
         gNMDA_deepng  (i,j)    = 0.d0
         gGABA_A_deepng  (i,j)  = 0.d0
          end do
          end do

         do L = firstcell, lastcell

c Handle LOT  -> deepng
      do i = 1, num_LOT_to_deepng   
       j = map_LOT_to_deepng(i,L) ! j = presynaptic cell
       k = com_LOT_to_deepng(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_LOT(j)  ! enumerate presyn. spikes
        presyntime = outtime_LOT(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_LOT_to_deepng  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepng(k,L)  = gAMPA_deepng(k,L) +
     &  gAMPA_LOT_to_deepng * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_LOT_to_deepng * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_LOT_to_deepng  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_LOT_to_deepng * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_LOT_to_deepng  
       if (gNMDA_deepng(k,L).gt.z)
     &  gNMDA_deepng(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle semilunar     -> deepng
      do i = 1, num_semilunar_to_deepng   
       j = map_semilunar_to_deepng(i,L) ! j = presynaptic cell
       k = com_semilunar_to_deepng(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_semilunar(j)  ! enumerate presyn. spikes
        presyntime = outtime_semilunar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_semilunar_to_deepng  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepng(k,L)  = gAMPA_deepng(k,L) +
     &  gAMPA_semilunar_to_deepng * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_semilunar_to_deepng * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_semilunar_to_deepng  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_semilunar_to_deepng * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_semilunar_to_deepng  
       if (gNMDA_deepng(k,L).gt.z)
     &  gNMDA_deepng(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle multipolar     -> deepng
      do i = 1, num_multipolar_to_deepng   
       j = map_multipolar_to_deepng(i,L) ! j = presynaptic cell
       k = com_multipolar_to_deepng(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_multipolar(j)  ! enumerate presyn. spikes
        presyntime = outtime_multipolar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_multipolar_to_deepng  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepng(k,L)  = gAMPA_deepng(k,L) +
     &  gAMPA_multipolar_to_deepng * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_multipolar_to_deepng * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_multipolar_to_deepng  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_multipolar_to_deepng * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_multipolar_to_deepng  
       if (gNMDA_deepng(k,L).gt.z)
     &  gNMDA_deepng(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask   -> deepng  
      do i = 1, num_deepbask_to_deepng      
       j = map_deepbask_to_deepng  (i,L) ! j = presynaptic cell
       k = com_deepbask_to_deepng  (i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_deepng      
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepng  (k,L)  = gGABA_A_deepng  (k,L) +
     &  gGABA_deepbask_to_deepng   * z      
! end GABA-A part

       end do ! m
      end do ! i

c Handle deepng     -> deepng
      do i = 1, num_deepng_to_deepng    
       j = map_deepng_to_deepng(i,L) ! j = presynaptic cell
       k = com_deepng_to_deepng(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepng(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepng(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepng_to_deepng    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepng(k,L)  = gGABA_A_deepng(k,L) +
     &  gGABA_deepng_to_deepng * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle placeholder5        -> deepng
      do i = 1, num_placeholder5_to_deepng 
       j = map_placeholder5_to_deepng(i,L) ! j = presynaptic cell
       k = com_placeholder5_to_deepng(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder5(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder5(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_placeholder5_to_deepng 
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepng(k,L)  = gAMPA_deepng(k,L) +
     &  gAMPA_placeholder5_to_deepng * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_placeholder5_to_deepng * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_placeholder5_to_deepng 
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_placeholder5_to_deepng * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_placeholder5_to_deepng  
       if (gNMDA_deepng(k,L).gt.z)
     &  gNMDA_deepng(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle L2pyr  -> deepng
      do i = 1, num_L2pyr_to_deepng
       j = map_L2pyr_to_deepng(i,L) ! j = presynaptic cell
       k = com_L2pyr_to_deepng(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L2pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L2pyr_to_deepng
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepng(k,L)  = gAMPA_deepng(k,L) +
     &  gAMPA_L2pyr_to_deepng * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_L2pyr_to_deepng * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L2pyr_to_deepng
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_L2pyr_to_deepng * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L2pyr_to_deepng
       if (gNMDA_deepng(k,L).gt.z)
     &  gNMDA_deepng  (k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle L3pyr  -> deepng
      do i = 1, num_L3pyr_to_deepng
       j = map_L3pyr_to_deepng(i,L) ! j = presynaptic cell
       k = com_L3pyr_to_deepng(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L3pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L3pyr_to_deepng
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepng(k,L)  = gAMPA_deepng(k,L) +
     &  gAMPA_L3pyr_to_deepng * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_L3pyr_to_deepng * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L3pyr_to_deepng
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepng(k,L) = gNMDA_deepng(k,L) +
     &  gNMDA_L3pyr_to_deepng * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L3pyr_to_deepng
       if (gNMDA_deepng(k,L).gt.z)
     &  gNMDA_deepng  (k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of deepng      
         ENDIF ! if (mod(O,how_often).eq.0) ...

! Define currents to deepng       cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for deepng     cells
       CALL INTEGRATE_deepng  (O, time, num_deepng ,
     &    V_deepng , curr_deepng ,
     & initialize, firstcell, lastcell,
     & gAMPA_deepng, gNMDA_deepng, gGABA_A_deepng,
     & Mg, 
     & gapcon_deepng  ,totSDgj_deepng   ,gjtable_deepng, dt,
     &  chi_deepng,mnaf_deepng,mnap_deepng,
     &  hnaf_deepng,mkdr_deepng,mka_deepng,
     &  hka_deepng,mk2_deepng,hk2_deepng,
     &  mkm_deepng,mkc_deepng,mkahp_deepng,
     &  mcat_deepng,hcat_deepng,mcal_deepng,
     &  mar_deepng)


        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
       do L = 1, num_deepng    
        distal_axon_deepintern (L + 300      ) = V_deepng     (59,L)
       end do
  
c          call mpi_allgather (distal_axon_deepng,
c    &  maxcellspernode, mpi_double_precision,
c    &  distal_axon_global,maxcellspernode,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)

        field_sup_local(1) = 0.d0     
        field_deep_local(1) = 0.d0     
c          call mpi_allgather (field_sup_local,     
c    &  1              , mpi_double_precision,
c    &  field_sup_global  , 1             ,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)
c          call mpi_allgather (field_deep_local,     
c    &  1              , mpi_double_precision,
c    &  field_deep_global  , 1             ,mpi_double_precision,
c    &                      MPI_COMM_WORLD, info)
  
           ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for deepng  

c      ELSE IF (THISNO.EQ.10) THEN
c      ELSE IF (nodecell(thisno) .eq. 'deepLTS ') THEN
c deepLTS

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 
          lastcell = num_deepLTS                            

        IF (mod(O,how_often).eq.0) then
c 1st set deepLTS  synaptic conductances to 0:

          do i = 1, numcomp_deepLTS
          do j = firstcell, lastcell
         gAMPA_deepLTS(i,j)    = 0.d0
         gNMDA_deepLTS(i,j)    = 0.d0
         gGABA_A_deepLTS(i,j)  = 0.d0 
          end do
          end do

         do L = firstcell, lastcell
c Handle L2pyr   -> deepLTS
      do i = 1, num_L2pyr_to_deepLTS  
       j = map_L2pyr_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_L2pyr_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L2pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L2pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L2pyr_to_deepLTS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepLTS(k,L)  = gAMPA_deepLTS(k,L) +
     &  gAMPA_L2pyr_to_deepLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_L2pyr_to_deepLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L2pyr_to_deepLTS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_L2pyr_to_deepLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L2pyr_to_deepLTS  
       if (gNMDA_deepLTS(k,L).gt.z)
     &  gNMDA_deepLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle placeholder3     -> deepLTS
      do i = 1, num_placeholder3_to_deepLTS    
       j = map_placeholder3_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_placeholder3_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder3(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder3(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_placeholder3_to_deepLTS    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepLTS(k,L)  = gGABA_A_deepLTS(k,L) +
     &  gGABA_placeholder3_to_deepLTS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle LOT  -> deepLTS
      do i = 1, num_LOT_to_deepLTS   
       j = map_LOT_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_LOT_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_LOT(j)  ! enumerate presyn. spikes
        presyntime = outtime_LOT(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_LOT_to_deepLTS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepLTS(k,L)  = gAMPA_deepLTS(k,L) +
     &  gAMPA_LOT_to_deepLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_LOT_to_deepLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_LOT_to_deepLTS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_LOT_to_deepLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_LOT_to_deepLTS  
       if (gNMDA_deepLTS(k,L).gt.z)
     &  gNMDA_deepLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle semilunar     -> deepLTS
      do i = 1, num_semilunar_to_deepLTS   
       j = map_semilunar_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_semilunar_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_semilunar(j)  ! enumerate presyn. spikes
        presyntime = outtime_semilunar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_semilunar_to_deepLTS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepLTS(k,L)  = gAMPA_deepLTS(k,L) +
     &  gAMPA_semilunar_to_deepLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_semilunar_to_deepLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_semilunar_to_deepLTS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_semilunar_to_deepLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_semilunar_to_deepLTS  
       if (gNMDA_deepLTS(k,L).gt.z)
     &  gNMDA_deepLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle multipolar     -> deepLTS
      do i = 1, num_multipolar_to_deepLTS   
       j = map_multipolar_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_multipolar_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_multipolar(j)  ! enumerate presyn. spikes
        presyntime = outtime_multipolar(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_multipolar_to_deepLTS  
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepLTS(k,L)  = gAMPA_deepLTS(k,L) +
     &  gAMPA_multipolar_to_deepLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_multipolar_to_deepLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_multipolar_to_deepLTS  
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_multipolar_to_deepLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_multipolar_to_deepLTS  
       if (gNMDA_deepLTS(k,L).gt.z)
     &  gNMDA_deepLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle deepbask   -> deepLTS
      do i = 1, num_deepbask_to_deepLTS    
       j = map_deepbask_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_deepbask_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_deepbask(j)  ! enumerate presyn. spikes
        presyntime = outtime_deepbask(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_deepbask_to_deepLTS    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepLTS(k,L)  = gGABA_A_deepLTS(k,L) +
     &  gGABA_deepbask_to_deepLTS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle supVIP    -> deepLTS
      do i = 1, num_supVIP_to_deepLTS    
       j = map_supVIP_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_supVIP_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_supVIP(j)  ! enumerate presyn. spikes
        presyntime = outtime_supVIP(m,j)
        delta = time - presyntime

! GABA-A part
        dexparg = delta / tauGABA_supVIP_to_deepLTS    
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gGABA_A_deepLTS(k,L)  = gGABA_A_deepLTS(k,L) +
     &  gGABA_supVIP_to_deepLTS * z      
! end GABA-A part

       end do ! m
      end do ! i


c Handle placeholder5        -> deepLTS
      do i = 1, num_placeholder5_to_deepLTS 
       j = map_placeholder5_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_placeholder5_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder5(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder5(m,j)
        delta = time - presyntime - thal_cort_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_placeholder5_to_deepLTS 
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepLTS(k,L)  = gAMPA_deepLTS(k,L) +
     &  gAMPA_placeholder5_to_deepLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_placeholder5_to_deepLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_placeholder5_to_deepLTS 
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_placeholder5_to_deepLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_placeholder5_to_deepLTS  
       if (gNMDA_deepLTS(k,L).gt.z)
     &  gNMDA_deepLTS(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


c Handle L3pyr  -> deepLTS
      do i = 1, num_L3pyr_to_deepLTS
       j = map_L3pyr_to_deepLTS(i,L) ! j = presynaptic cell
       k = com_L3pyr_to_deepLTS(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L3pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_L3pyr_to_deepLTS
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_deepLTS(k,L)  = gAMPA_deepLTS(k,L) +
     &  gAMPA_L3pyr_to_deepLTS * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_L3pyr_to_deepLTS * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L3pyr_to_deepLTS
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_deepLTS(k,L) = gNMDA_deepLTS(k,L) +
     &  gNMDA_L3pyr_to_deepLTS * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L3pyr_to_deepLTS
       if (gNMDA_deepLTS(k,L).gt.z)
     &  gNMDA_deepLTS(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


         end do
c End enumeration of deepLTS    
        ENDIF  !  if (mod(O,how_often).eq.0) ...

! Define currents to deepLTS     cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for deepLTS     cells
       CALL INTEGRATE_deepLTS (O, time, num_deepLTS ,
     &    V_deepLTS , curr_deepLTS ,
     & initialize, firstcell, lastcell,
     & gAMPA_deepLTS, gNMDA_deepLTS, gGABA_A_deepLTS,
     & Mg, 
     & gapcon_deepLTS  ,totSDgj_deepLTS   ,gjtable_deepLTS, dt,
     &  chi_deepLTS,mnaf_deepLTS,mnap_deepLTS,
     &  hnaf_deepLTS,mkdr_deepLTS,mka_deepLTS,
     &  hka_deepLTS,mk2_deepLTS,hk2_deepLTS,
     &  mkm_deepLTS,mkc_deepLTS,mkahp_deepLTS,
     &  mcat_deepLTS,hcat_deepLTS,mcal_deepLTS,
     &  mar_deepLTS)


        IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
       do L = 1, num_deepLTS    
        distal_axon_deepintern   (L + 200      ) = V_deepLTS     (59,L)
       end do
  
           call mpi_allgather (distal_axon_deepintern,
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_sup_local(1) = 0.d0     
        field_deep_local(1) = 0.d0     
           call mpi_allgather (field_sup_local,     
     &  1              , mpi_double_precision,
     &  field_sup_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_deep_local,     
     &  1              , mpi_double_precision,
     &  field_deep_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
        ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for deepLTS


c      ELSE IF (THISNO.EQ.12) THEN
       ELSE IF (nodecell(thisno) .eq. 'placeholder5') THEN
c placeholder5

c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 
          lastcell = num_placeholder5                               

        IF (mod(O,how_often).eq.0) then
c 1st set placeholder5 synaptic conductances to 0:

          do i = 1, numcomp_placeholder5
          do j = firstcell, lastcell
         gAMPA_placeholder5(i,j)         = 0.d0 
         gNMDA_placeholder5(i,j)         = 0.d0
         gGABA_A_placeholder5(i,j)       = 0.d0 
         gGABA_B_placeholder5(i,j)       = 0.d0 
          end do
          end do

         do L = firstcell, lastcell
c Handle placeholder6       -> placeholder5
      do i = 1, num_placeholder6_to_placeholder5     
       j = map_placeholder6_to_placeholder5(i,L) ! j = presynaptic cell
       k = com_placeholder6_to_placeholder5(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder6(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder6(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part
        dexparg1 = delta / tauGABA1_placeholder6_to_placeholder5     
c note that dexparg1 = MINUS the actual arg. to dexp
         if (dexparg1.le.5.d0) then
          z1 = dexptablesmall (int(dexparg1*1000.d0))
         else if (dexparg1.le.100.d0) then
          z1 = dexptablebig (int(dexparg1*10.d0))
         else
          z1 = 0.d0
         endif

        dexparg2 = delta / tauGABA2_placeholder6_to_placeholder5     
c note that dexparg2 = MINUS the actual arg. to dexp
         if (dexparg2.le.5.d0) then
          z2 = dexptablesmall (int(dexparg2*1000.d0))
         else if (dexparg2.le.100.d0) then
          z2 = dexptablebig (int(dexparg2*10.d0))
         else
          z2 = 0.d0
         endif

      gGABA_A_placeholder5(k,L)  = gGABA_A_placeholder5(k,L) +
     &  gGABA_placeholder6_to_placeholder5*(0.625d0*z1+0.375d0*z2) 
! end GABA-A part


      gGABA_B_placeholder5(k,L) = gGABA_B_placeholder5(k,L) +
     &   gGABAB_placeholder6_to_placeholder5 * otis_table(k0)
! end GABA-B part
       end do ! m
      end do ! i


c Handle L3pyr -> placeholder5
      do i = 1, num_L3pyr_to_placeholder5
       j = map_L3pyr_to_placeholder5(i,L) ! j = presynaptic cell
       k = com_L3pyr_to_placeholder5(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L3pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime - cort_thal_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_L3pyr_to_placeholder5
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder5(k,L)  = gAMPA_placeholder5(k,L) +
     &  gAMPA_L3pyr_to_placeholder5 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder5(k,L) = gNMDA_placeholder5(k,L) +
     &  gNMDA_L3pyr_to_placeholder5 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L3pyr_to_placeholder5
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder5(k,L) = gNMDA_placeholder5(k,L) +
     &  gNMDA_L3pyr_to_placeholder5 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L3pyr_to_placeholder5 
       if (gNMDA_placeholder5(k,L).gt.z)
     &  gNMDA_placeholder5(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


         end do
c End enumeration of placeholder5         
          ENDIF  !  if (mod(O,how_often).eq.0) ...

! Define currents to placeholder5          cells, ectopic spikes,
! tonic synaptic conductances

      if (mod(O,200).eq.0) then
       call durand(seed,num_placeholder5     ,ranvec_placeholder5     ) 
        do L = firstcell, lastcell
         if ((ranvec_placeholder5     (L).gt.0.d0).and.
     &  (ranvec_placeholder5     (L).le.noisepe_placeholder5     )) then
          curr_placeholder5     (135,L) = 0.4d0
         else
          curr_placeholder5     (135,L) = 0.d0
         endif 
        end do
      endif

! Call integration routine for placeholder5          cells
       CALL INTEGRATE_placeholder5     (O, time, num_placeholder5      ,
     &    V_placeholder5      , curr_placeholder5      ,
     & initialize, firstcell, lastcell,
     & gAMPA_placeholder5, gNMDA_placeholder5, gGABA_A_placeholder5 ,
     & gGABA_B_placeholder5, Mg, 
     & gapcon_placeholder5,totaxgj_placeholder5,gjtable_placeholder5,dt,
     &  chi_placeholder5,mnaf_placeholder5,mnap_placeholder5,
     &  hnaf_placeholder5,mkdr_placeholder5,mka_placeholder5,
     &  hka_placeholder5,mk2_placeholder5,hk2_placeholder5,
     &  mkm_placeholder5,mkc_placeholder5,mkahp_placeholder5,
     &  mcat_placeholder5,hcat_placeholder5,mcal_placeholder5,
     &  mar_placeholder5)
9144    CONTINUE


         IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_placeholder5         
       do L = firstcell, lastcell
        distal_axon_placeholder5 (L-firstcell+1) = V_placeholder5(135,L)
       end do
  
           call mpi_allgather (distal_axon_placeholder5,     
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_sup_local(1) = 0.d0       
        field_deep_local(1) = 0.d0     
           call mpi_allgather (field_sup_local,     
     &  1              , mpi_double_precision,
     &  field_sup_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_deep_local,     
     &  1              , mpi_double_precision,
     &  field_deep_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
        ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for placeholder5

c      ELSE IF (THISNO.EQ.13) THEN
       ELSE IF (nodecell(thisno) .eq. 'placeholder6') THEN
c placeholder6

c Determine which particular cells this node will be concerned with.
c         i = place (thisno)
          firstcell = 1 
          lastcell = num_placeholder6                               

        IF (mod(O,how_often).eq.0) then
c 1st set placeholder6 synaptic conductances to 0:

          do i = 1, numcomp_placeholder6
          do j = firstcell, lastcell
         gAMPA_placeholder6(i,j)         = 0.d0 
         gNMDA_placeholder6(i,j)         = 0.d0
         gGABA_A_placeholder6(i,j)       = 0.d0
         gGABA_B_placeholder6(i,j)       = 0.d0
          end do
          end do

         do L = firstcell, lastcell
c Handle placeholder5        -> placeholder6
      do i = 1, num_placeholder5_to_placeholder6
       j = map_placeholder5_to_placeholder6(i,L) ! j = presynaptic cell
       k = com_placeholder5_to_placeholder6(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder5(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder5(m,j)
        delta = time - presyntime

! AMPA part
        dexparg = delta / tauAMPA_placeholder5_to_placeholder6
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder6(k,L)  = gAMPA_placeholder6(k,L) +
     &  gAMPA_placeholder5_to_placeholder6 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder6(k,L) = gNMDA_placeholder6(k,L) +
     &  gNMDA_placeholder5_to_placeholder6 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_placeholder5_to_placeholder6 
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder6(k,L) = gNMDA_placeholder6(k,L) +
     &  gNMDA_placeholder5_to_placeholder6 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_placeholder5_to_placeholder6
       if (gNMDA_placeholder6(k,L).gt.z)
     &  gNMDA_placeholder6(k,L) = z
! end NMDA part

       end do ! m
      end do ! i


c Handle placeholder6        -> placeholder6
      do i = 1, num_placeholder6_to_placeholder6     
       j = map_placeholder6_to_placeholder6(i,L) ! j = presynaptic cell
       k = com_placeholder6_to_placeholder6(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_placeholder6(j)  ! enumerate presyn. spikes
        presyntime = outtime_placeholder6(m,j)
        delta = time - presyntime
        k0 = nint (10.d0 * delta) ! time, in units of 0.1 ms, to pass to otis_table
        if (k0 .gt. 50000) k = 50000  ! limit on size of otis_table

! GABA-A part
        dexparg1 = delta / tauGABA1_placeholder6_to_placeholder6     
c note that dexparg1 = MINUS the actual arg. to dexp
         if (dexparg1.le.5.d0) then
          z1 = dexptablesmall (int(dexparg1*1000.d0))
         else if (dexparg1.le.100.d0) then
          z1 = dexptablebig (int(dexparg1*10.d0))
         else
          z1 = 0.d0
         endif

        dexparg2 = delta / tauGABA2_placeholder6_to_placeholder6     
c note that dexparg2 = MINUS the actual arg. to dexp
         if (dexparg2.le.5.d0) then
          z2 = dexptablesmall (int(dexparg2*1000.d0))
         else if (dexparg2.le.100.d0) then
          z2 = dexptablebig (int(dexparg2*10.d0))
         else
          z2 = 0.d0
         endif

      gGABA_A_placeholder6(k,L)  = gGABA_A_placeholder6(k,L) +
     & gGABA_placeholder6_to_placeholder6 * (0.56d0 * z1 + 0.44d0 * z2) 
! end GABA-A part

      gGABA_B_placeholder6(k,L) = gGABA_B_placeholder6(k,L) +
     &   gGABAB_placeholder6_to_placeholder6 * otis_table(k0)

! end GABA-B part
       end do ! m
      end do ! i


c Handle L3pyr  -> placeholder6
      do i = 1, num_L3pyr_to_placeholder6
       j = map_L3pyr_to_placeholder6(i,L) ! j = presynaptic cell
       k = com_L3pyr_to_placeholder6(i,L) ! k = comp. on postsyn. cell

       do m = 1, outctr_L3pyr(j)  ! enumerate presyn. spikes
        presyntime = outtime_L3pyr(m,j)
        delta = time - presyntime - cort_thal_delay

         IF (DELTA.GE.0.d0) THEN
! AMPA part
        dexparg = delta / tauAMPA_L3pyr_to_placeholder6
c note that dexparg = MINUS the actual arg. to dexp
         if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif

      gAMPA_placeholder6(k,L)  = gAMPA_placeholder6(k,L) +
     &  gAMPA_L3pyr_to_placeholder6 * delta * z      
! end AMPA part

! NMDA part
        if (delta.le.5.d0) then
       gNMDA_placeholder6(k,L) = gNMDA_placeholder6(k,L) +
     &  gNMDA_L3pyr_to_placeholder6 * delta * 0.2d0
        else
       dexparg = (delta - 5.d0)/tauNMDA_L3pyr_to_placeholder6
          if (dexparg.le.5.d0) then
          z = dexptablesmall (int(dexparg*1000.d0))
         else if (dexparg.le.100.d0) then
          z = dexptablebig (int(dexparg*10.d0))
         else
          z = 0.d0
         endif
       gNMDA_placeholder6(k,L) = gNMDA_placeholder6(k,L) +
     &  gNMDA_L3pyr_to_placeholder6 * z
        endif
c Test for NMDA saturation
       z = NMDA_saturation_fact * gNMDA_L3pyr_to_placeholder6 
       if (gNMDA_placeholder6(k,L).gt.z)
     &  gNMDA_placeholder6(k,L) = z
! end NMDA part

        ENDIF  ! condition for checking that delta >= 0.
       end do ! m
      end do ! i


         end do
c End enumeration of placeholder6         
        ENDIF  !  if (mod(O,how_often).eq.0) ...

! Define currents to placeholder6          cells, ectopic spikes,
! tonic synaptic conductances

! Call integration routine for placeholder6          cells
       CALL INTEGRATE_placeholder6  (O, time, num_placeholder6      ,
     &    V_placeholder6      , curr_placeholder6      ,
     & initialize, firstcell, lastcell,
     & gAMPA_placeholder6, gNMDA_placeholder6, gGABA_A_placeholder6  ,
     & gGABA_B_placeholder6, Mg, 
     & gapcon_placeholder6,totaxgj_placeholder6,gjtable_placeholder6,dt,
     &  chi_placeholder6,mnaf_placeholder6,mnap_placeholder6,
     &  hnaf_placeholder6,mkdr_placeholder6,mka_placeholder6,
     &  hka_placeholder6,mk2_placeholder6,hk2_placeholder6,
     &  mkm_placeholder6,mkc_placeholder6,mkahp_placeholder6,
     &  mcat_placeholder6,hcat_placeholder6,mcal_placeholder6,
     &  mar_placeholder6,field_sup,field_deep,rel_axonshift_L2pyr)


         IF (mod(O,how_often).eq.0) then
! Set up distal axon voltage array and broadcast it.
c      do L = 1, num_placeholder6         
       do L = firstcell, lastcell
        distal_axon_placeholder6 (L-firstcell+1) = V_placeholder6(74,L)
       end do
  
           call mpi_allgather (distal_axon_placeholder6,      
     &  maxcellspernode, mpi_double_precision,
     &  distal_axon_global,maxcellspernode,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)

        field_sup_local(1) = 0.d0      
        field_deep_local(1) = 0.d0       
           call mpi_allgather (field_sup_local,     
     &  1              , mpi_double_precision,
     &  field_sup_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
           call mpi_allgather (field_deep_local,     
     &  1              , mpi_double_precision,
     &  field_deep_global  , 1             ,mpi_double_precision,
     &                      MPI_COMM_WORLD, info)
  
         ENDIF  !  if (mod(O,how_often).eq.0) ...

! END thisno for placeholder6

       ENDIF  ! if (mod(O,how_often).eq.0) then ...

! Update distal axon vectors, then outctr's and outtime tables.
! This code is common to all the nodes.
! Some of this section adapted from supergj.f
c     IF (mod(O,how_often).eq.0) then
      IF (mod(O,  5      ).eq.0) then ! Necessary because gj data also
!  being updated, not just synaptic
c Construct distal axon vectors, taking into account the structure of
c distal_axon_global: let m = maxcellspernode;
c then nodesfor_L2pyr segments, each m entries long;
! Do the same for voltages at sites of possible axonal gj - now obsolete

            ictr = 0 ! will keep track of which segment in distal_axon_global

c Make the unpacking "explicit"
            do L = 1, 1000 
              ldistal_axon_L2pyr(L) = distal_axon_global (L)
            end do
            do L = 1, 100
         ldistal_axon_placeholder1(L) = distal_axon_global (1000+L)
            end do
            do L = 1, 100
         ldistal_axon_placeholder2(L) = distal_axon_global (1100+L)
            end do
            do L = 1, 100
         ldistal_axon_placeholder3(L) = distal_axon_global (1200+L)
            end do
            do L = 1, 100
         ldistal_axon_supVIP  (L)  = distal_axon_global (1300+L)
            end do
            do L = 1, 100
         ldistal_axon_supng  (L)   = distal_axon_global (1400+L)
            end do
            do L = 1, 500
         ldistal_axon_LOT(L) = distal_axon_global (1500+L)
            end do
            do L = 1, 500
         ldistal_axon_semilunar (L)   = distal_axon_global (2000+L)
            end do
            do L = 1, 500
         ldistal_axon_multipolar  (L) = distal_axon_global (2500+L)
            end do 
            do L = 1, 500
         ldistal_axon_L3pyr(L) = distal_axon_global (3000+L)
            end do
            do L = 1, 200
         ldistal_axon_deepbask(L)  = distal_axon_global (3500+L)
            end do
            do L = 1, 100
         ldistal_axon_deepLTS(L) =  distal_axon_global (3700+L)
            end do
            do L = 1, 200
         ldistal_axon_deepng (L)  =   distal_axon_global (3800+L)
            end do
            do L = 1, 500
        ldistal_axon_placeholder5  (L)    =  distal_axon_global (4000+L)
            end do
            do L = 1, 500
         ldistal_axon_placeholder6  (L)   =  distal_axon_global (4500+L)
            end do

c End updating of distal axon vectors.

       do L = 1, num_L2pyr
        if (ldistal_axon_L2pyr(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
        if (outctr_L2pyr(L).eq.0) then
            outctr_L2pyr(L) = 1
            outtime_L2pyr(1,L) = time
          else
      if ((time-outtime_L2pyr(outctr_L2pyr(L),L))
     &   .gt. axon_refrac_time) then
           outctr_L2pyr(L) = outctr_L2pyr(L) + 1
           outtime_L2pyr (outctr_L2pyr(L),L) = time
            endif
          endif
       endif
       end do  ! do L = 1, num_L2pyr


       do L = 1, num_placeholder1   
        if (ldistal_axon_placeholder1(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
        if (outctr_placeholder1(L).eq.0) then
            outctr_placeholder1(L) = 1
            outtime_placeholder1(1,L) = time
          else
      if ((time-outtime_placeholder1(outctr_placeholder1(L),L))
     &   .gt. axon_refrac_time) then
             outctr_placeholder1(L) = outctr_placeholder1(L) + 1
            outtime_placeholder1 (outctr_placeholder1(L),L) = time
            endif
          endif
         endif
       end do  ! do L = 1, num_placeholder1  

       do L = 1, num_supng   
         if (ldistal_axon_supng(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
        if (outctr_supng(L).eq.0) then
         outctr_supng(L) = 1
         outtime_supng(1,L) = time
          else
      if ((time-outtime_supng(outctr_supng(L),L))
     &   .gt. axon_refrac_time) then
           outctr_supng(L) = outctr_supng(L) + 1
           outtime_supng (outctr_supng(L),L) = time
            endif
          endif
         endif
       end do  ! do L = 1, num_supng  

       do L = 1, num_placeholder2   
         if (ldistal_axon_placeholder2(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
         if (outctr_placeholder2(L).eq.0) then
             outctr_placeholder2(L) = 1
             outtime_placeholder2(1,L) = time
          else
      if ((time-outtime_placeholder2(outctr_placeholder2(L),L))
     &   .gt. axon_refrac_time) then
         outctr_placeholder2(L) = outctr_placeholder2(L) + 1
         outtime_placeholder2 (outctr_placeholder2(L),L) = time
            endif
          endif
        endif
       end do  ! do L = 1, num_placeholder2  

       do L = 1, num_placeholder3    
        if (ldistal_axon_placeholder3(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
        if (outctr_placeholder3(L).eq.0) then
            outctr_placeholder3(L) = 1
            outtime_placeholder3(1,L) = time
          else
      if ((time-outtime_placeholder3(outctr_placeholder3(L),L))
     &   .gt. axon_refrac_time) then
             outctr_placeholder3(L) = outctr_placeholder3(L) + 1
             outtime_placeholder3 (outctr_placeholder3(L),L) = time
            endif
          endif
          endif
       end do  ! do L = 1, num_placeholder3  

       do L = 1, num_LOT 
        if (ldistal_axon_LOT(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
         if (outctr_LOT(L).eq.0) then
            outctr_LOT(L) = 1
            outtime_LOT(1,L) = time
          else
      if ((time-outtime_LOT(outctr_LOT(L),L))
     &   .gt. axon_refrac_time) then
         outctr_LOT(L) = outctr_LOT(L) + 1
         outtime_LOT (outctr_LOT(L),L) = time
            endif
          endif
       endif
       end do  ! do L = 1, num_LOT

       do L = 1, num_semilunar    
        if (ldistal_axon_semilunar(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
        if (outctr_semilunar(L).eq.0) then
            outctr_semilunar(L) = 1
            outtime_semilunar(1,L) = time
          else
      if ((time-outtime_semilunar(outctr_semilunar(L),L))
     &   .gt. axon_refrac_time) then
         outctr_semilunar(L) = outctr_semilunar(L) + 1
         outtime_semilunar (outctr_semilunar(L),L) = time
            endif
          endif
       endif
       end do  ! do L = 1, num_semilunar   

       do L = 1, num_multipolar    
           if (ldistal_axon_multipolar(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
           if (outctr_multipolar(L).eq.0) then
            outctr_multipolar(L) = 1
            outtime_multipolar(1,L) = time
          else
      if ((time-outtime_multipolar(outctr_multipolar(L),L))
     &   .gt. axon_refrac_time) then
         outctr_multipolar(L) = outctr_multipolar(L) + 1
         outtime_multipolar (outctr_multipolar(L),L) = time
            endif
          endif
        endif
       end do  ! do L = 1, num_multipolar   

       do L = 1, num_L3pyr    
        if (ldistal_axon_L3pyr(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
          if (outctr_L3pyr(L).eq.0) then
            outctr_L3pyr(L) = 1
            outtime_L3pyr(1,L) = time
          else
      if ((time-outtime_L3pyr(outctr_L3pyr(L),L))
     &   .gt. axon_refrac_time) then
             outctr_L3pyr(L) = outctr_L3pyr(L) + 1
             outtime_L3pyr (outctr_L3pyr(L),L) = time
            endif
          endif
       endif
       end do  ! do L = 1, num_L3pyr   

       do L = 1, num_deepbask     
         if (ldistal_axon_deepbask(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
          if (outctr_deepbask(L).eq.0) then
            outctr_deepbask(L) = 1
            outtime_deepbask(1,L) = time
          else
      if ((time-outtime_deepbask(outctr_deepbask(L),L))
     &   .gt. axon_refrac_time) then
        outctr_deepbask(L) = outctr_deepbask(L) + 1
        outtime_deepbask (outctr_deepbask(L),L) = time
            endif
          endif
        endif
       end do  ! do L = 1, num_deepbask   

       do L = 1, num_deepng     
        if (ldistal_axon_deepng(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
          if (outctr_deepng(L).eq.0) then
            outctr_deepng(L) = 1
            outtime_deepng(1,L) = time
          else
      if ((time-outtime_deepng(outctr_deepng(L),L))
     &   .gt. axon_refrac_time) then
           outctr_deepng(L) = outctr_deepng(L) + 1
           outtime_deepng (outctr_deepng(L),L) = time
            endif
          endif
        endif
       end do  ! do L = 1, num_deepng   

       do L = 1, num_deepLTS     
         if (ldistal_axon_deepLTS(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
          if (outctr_deepLTS(L).eq.0) then
            outctr_deepLTS(L) = 1
            outtime_deepLTS(1,L) = time
          else
      if ((time-outtime_deepLTS(outctr_deepLTS(L),L))
     &   .gt. axon_refrac_time) then
         outctr_deepLTS(L) = outctr_deepLTS(L) + 1
          outtime_deepLTS (outctr_deepLTS(L),L) = time
            endif
          endif
        endif
       end do  ! do L = 1, num_deepLTS   

       do L = 1, num_supVIP      
          if (ldistal_axon_supVIP(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
          if (outctr_supVIP(L).eq.0) then
            outctr_supVIP(L) = 1
            outtime_supVIP(1,L) = time
          else
      if ((time-outtime_supVIP(outctr_supVIP(L),L))
     &   .gt. axon_refrac_time) then
             outctr_supVIP(L) = outctr_supVIP(L) + 1
             outtime_supVIP (outctr_supVIP(L),L) = time
            endif
          endif
       endif
       end do  ! do L = 1, num_supVIP   

       do L = 1, num_placeholder5      
        if (ldistal_axon_placeholder5(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
         if (outctr_placeholder5(L).eq.0) then
            outctr_placeholder5(L) = 1
            outtime_placeholder5(1,L) = time
          else
      if ((time-outtime_placeholder5(outctr_placeholder5(L),L))
     &   .gt. axon_refrac_time) then
         outctr_placeholder5(L) = outctr_placeholder5(L) + 1
         outtime_placeholder5 (outctr_placeholder5(L),L) = time
            endif
          endif
        endif
       end do  ! do L = 1, num_placeholder5   

       do L = 1, num_placeholder6      
        if (ldistal_axon_placeholder6(L).ge.0.d0) then
c with threshold = 0, means axonal spike must be overshooting.
         if (outctr_placeholder6(L).eq.0) then
            outctr_placeholder6(L) = 1
            outtime_placeholder6(1,L) = time
          else
      if ((time-outtime_placeholder6(outctr_placeholder6(L),L))
     &   .gt. axon_refrac_time) then
             outctr_placeholder6(L) = outctr_placeholder6(L) + 1
           outtime_placeholder6 (outctr_placeholder6(L),L) = time
            endif
          endif
        endif
       end do  ! do L = 1, num_placeholder6   

       field_sup_tot = 0.d0
       field_deep_tot = 0.d0
        do i = 1, numnodes
         field_sup_tot = field_sup_tot + field_sup_global(i)
         field_deep_tot = field_deep_tot + field_deep_global(i)
        end do

      ENDIF  ! if (mod(O,how_often).eq.0) ...
       ! CHANGED to if (mod(O,5).eq.0)...
! End updating outctr's and outtime tables, and computing fields


! Set up output data to be written
c        GOTO 2000 ! for testing
       if (mod(O, 50) == 0) then


c      if (thisno.eq.0) then
       IF (nodecell(thisno) .eq. 'L2pyr       ') THEN
c L2pyr
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
           if (i.eq.1) then
          firstcell = 1 
           else
          firstcell = 501
           endif
          lastcell = firstcell -1 +  ncellspernode

        outrcd( 1) = time
        outrcd( 2) = v_L2pyr(1,firstcell+1)
        outrcd( 3) = v_L2pyr(numcomp_L2pyr,firstcell+1)
        outrcd( 4) = v_L2pyr(43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_L2pyr(1,i)
          end do
        outrcd( 5) = z / dble(lastcell - firstcell + 1) ! - av. cell somata 
         z = 0.d0
          do i = 1, numcomp_L2pyr
           z = z + gAMPA_L2pyr(i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_L2pyr
           z = z + gNMDA_L2pyr(i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_L2pyr
           z = z + gGABA_A_L2pyr(i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
         z = 0.d0
          do i = 1, numcomp_L2pyr
           z = z + gGABA_B_L2pyr(i,firstcell+1)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 2 
        outrcd(10) = v_L2pyr(1, 426       )
        outrcd(11) = v_L2pyr(1,firstcell+3)
          z = 0.d0
          do i = firstcell, lastcell
           if(v_L2pyr(numcomp_L2pyr,i) .gt. 0.d0) z = z + 1.d0
          end do
        outrcd(12) = z   
        outrcd(13) = 0.d0 ! field_sup_tot     
        outrcd(14) = 0.d0 ! field_deep_tot     
        outrcd(15) = v_L2pyr(1, firstcell+21     )
        outrcd(16) = v_L2pyr(43,firstcell+21     )
        outrcd(17) = v_L2pyr(1,firstcell+387)
        outrcd(18) = v_L2pyr(43,firstcell+387)

            if (place(thisno).eq.1) then
      OPEN(11,FILE='piriformENDO96.L2pyr')
      WRITE (11,FMT='(18F10.4)') (OUTRCD(I),I=1,18)
c           else
         outrcd( 1) = time
         outrcd( 2) = V_L2pyr( 1, 416)
         outrcd( 3) = V_L2pyr( 1, 160)
         outrcd( 4) = V_L2pyr( 1, 310)
         outrcd( 5) = V_L2pyr( 1, 343)
         outrcd( 6) = V_L2pyr( 1,  93)
         outrcd( 7) = V_L2pyr( 1,  65)
         outrcd( 8) = V_L2pyr( 1,  86)
         outrcd( 9) = V_L2pyr( 1, 154)
         outrcd(10) = V_L2pyr( 1, 466)
         outrcd(11) = V_L2pyr( 1,  93)
         outrcd(12) = V_L2pyr( 1, 102)
         outrcd(13) = V_L2pyr( 1, 158)
         outrcd(14) = V_L2pyr( 1,  22)
         outrcd(15) = V_L2pyr( 1,  40)
         outrcd(16) = V_L2pyr( 1,  60)
         outrcd(17) = V_L2pyr( 1,  66)
         outrcd(18) = V_L2pyr( 1, 216)
         outrcd(19) = V_L2pyr( 1, 320)
         outrcd(20) = V_L2pyr( 1, 328)
      OPEN(111,FILE='piriformENDO96.L2pyrA')
      WRITE (111,FMT='(20F10.4)') (OUTRCD(I),I=1,20)
            end if

       do L = firstcell, lastcell     
c       if (v_L2pyr (1,L) .ge. -15.d0) then
        if (v_L2pyr (1,L) .ge. -25.d0) then
          if (place(thisno).eq.1) then
         OPEN(41,FILE='piriformENDO96.L2pyrrast')
         WRITE (41,8789) time, L
          else
         OPEN(411,FILE='piriformENDO96.L2pyrrastA')
         WRITE (411,8789) time, L
          endif
8789     FORMAT (f8.2,3x,i5)
        end if
        if (v_L2pyr (numcomp_L2pyr,L) .ge. 0.d0) then
          if (place(thisno).eq.1) then
         OPEN(42,FILE='piriformENDO96.L2pyrrastax')
         WRITE (42,8789) time, L
          else
         OPEN(421,FILE='piriformENDO96.L2pyrrastaxA')
         WRITE (421,8789) time, L
          endif
c This only records the 1st 500 L2pyr cells?
        end if
       end do

       else if (thisno.eq.2) then
c      else IF (nodecell(thisno) .eq. 'supintern   ') THEN
c placeholder1 
c Determine which particular cells this node will be concerned with.
          i = place (thisno)
          firstcell = 1 
          lastcell = num_placeholder1                           

        outrcd( 1) = time
        outrcd( 2) = v_placeholder1  (1,firstcell+1)
        outrcd( 3)=v_placeholder1  (numcomp_placeholder1,firstcell+1)
        outrcd( 4) = v_placeholder1  (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_placeholder1(1,i)
          end do
        outrcd( 5) = z / dble(lastcell - firstcell + 1  ) ! - av. cell somata 
         z = 0.d0
          do i = 1, numcomp_placeholder1   
           z = z + gAMPA_placeholder1  (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_placeholder1   
           z = z + gNMDA_placeholder1  (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_placeholder1  
           z = z + gGABA_A_placeholder1  (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_placeholder1  (1,firstcell+2)
        outrcd(10) = v_placeholder1  (1,firstcell+3)
c     OPEN(13,FILE='piriformENDO96.placeholder1')
c     WRITE (13,FMT='(10F10.4)') (OUTRCD(I),I=1,10)

c supng 
c Determine which particular cells this node will be concerned with.
          firstcell = 1 
          lastcell = num_supng                           

        outrcd( 1) = time
        outrcd( 2) = v_supng  (1,firstcell+1)
        outrcd( 3) = v_supng  (numcomp_supng,firstcell+1)
        outrcd( 4) = v_supng  (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_supng(1,i)
          end do
        outrcd( 5) = z / dble(lastcell - firstcell + 1  ) ! - av. cell somata 
         z = 0.d0
          do i = 1, numcomp_supng   
           z = z + gAMPA_supng  (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_supng   
           z = z + gNMDA_supng  (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_supng  
           z = z + gGABA_A_supng  (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_supng  (1,firstcell+2)
        outrcd(10) = v_supng  (1,firstcell+3)
         if (place(thisno).eq.1) then
      OPEN(33,FILE='piriformENDO96.supng')
      WRITE (33,FMT='(10F10.4)') (OUTRCD(I),I=1,10)
         endif

       do L = firstcell, lastcell     
        if (v_supng     (1,L) .ge. -25.d0) then
          if (place(thisno).eq.1) then
         OPEN(816,FILE='piriformENDO96.supngrast   ')
         WRITE (816,8789) time, L
          endif
        endif
       end do

c placeholder2 
c Determine which particular cells this node will be concerned with.
          firstcell = 1 
          lastcell = num_placeholder2                              

        outrcd( 1) = time
        outrcd( 2) = v_placeholder2 (1,firstcell+1)
        outrcd( 3) = v_placeholder2 (numcomp_placeholder2 ,firstcell+1)
        outrcd( 4) = v_placeholder2  (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_placeholder2(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_placeholder2  
           z = z + gAMPA_placeholder2  (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_placeholder2  
           z = z + gNMDA_placeholder2  (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_placeholder2  
           z = z + gGABA_A_placeholder2  (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_placeholder2  (1,firstcell+2)
        outrcd(10) = v_placeholder2  (1,firstcell+3)
          if (place(thisno).eq.1) then
c     OPEN(14,FILE='piriformENDO96.placeholder2')
c     WRITE (14,FMT='(10F10.4)') (OUTRCD(I),I=1,10)
          endif

c placeholder3  
c Determine which particular cells this node will be concerned with.
          firstcell = 1 
          lastcell = num_placeholder3                            

        outrcd( 1) = time
        outrcd( 2) = v_placeholder3   (1,firstcell+1)
        outrcd( 3) = v_placeholder3   (numcomp_placeholder3,firstcell+1)
        outrcd( 4) = v_placeholder3   (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_placeholder3(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_placeholder3   
           z = z + gAMPA_placeholder3   (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_placeholder3   
           z = z + gNMDA_placeholder3   (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_placeholder3   
           z = z + gGABA_A_placeholder3   (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_placeholder3   (1,firstcell+2)
        outrcd(10) = v_placeholder3   (1,firstcell+3)
c     OPEN(15,FILE='piriformENDO96.placeholder3')
c     WRITE (15,FMT='(10F10.4)') (OUTRCD(I),I=1,10)

c supVIP 
c Determine which particular cells this node will be concerned with.
          firstcell = 1 
          lastcell = num_supVIP 

        outrcd( 1) = time
        outrcd( 2) = v_supVIP  (1,firstcell+1)
        outrcd( 3) = v_supVIP  (numcomp_supVIP  ,firstcell+1)
        outrcd( 4) = v_supVIP  (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_supVIP(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_supVIP  
           z = z + gAMPA_supVIP  (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_supVIP   
           z = z + gNMDA_supVIP  (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_supVIP  
           z = z + gGABA_A_supVIP  (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_supVIP  (1,firstcell+2)
        outrcd(10) = v_supVIP  (1,firstcell+3)
      OPEN(22,FILE='piriformENDO96.supVIP')
      WRITE (22,FMT='(10F10.4)') (OUTRCD(I),I=1,10)


       do L = firstcell, lastcell     
        if (v_supVIP    (1,L) .ge. -25.d0) then
          if (place(thisno).eq.1) then
         OPEN(616,FILE='piriformENDO96.supVIPrast   ')
         WRITE (616,8789) time, L
          endif
        endif
       end do

       else IF (nodecell(thisno) .eq. 'LOT         ') THEN
c LOT
c Determine which particular cells this node will be concerned with.
          firstcell = 1 
          lastcell = num_LOT 

        outrcd( 1) = time
        outrcd( 2) = v_LOT(1,firstcell+1)
        outrcd( 3) = v_LOT(numcomp_LOT,firstcell+1)
        outrcd( 4) = v_LOT(43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_LOT(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_LOT
           z = z + gAMPA_LOT(i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 1 
         z = 0.d0
          do i = 1, numcomp_LOT
           z = z + gNMDA_LOT(i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 1 
         z = 0.d0
          do i = 1, numcomp_LOT
           z = z + gGABA_A_LOT(i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 1 
         z = 0.d0
          do i = 1, numcomp_LOT
           z = z + gGABA_B_LOT(i,firstcell+1)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 1 
        outrcd(10) = v_LOT(1, 201       )
        outrcd(11) = v_LOT(1, 250       )

          if (place(thisno).eq.1) then
      OPEN(16,FILE='piriformENDO96.LOT')
      WRITE (16,FMT='(11F10.4)') (OUTRCD(I),I=1,11)
          endif

       do L = firstcell, lastcell     
        if (v_LOT(1,L) .ge. -25.d0) then
          if (place(thisno).eq.1) then
         OPEN(516,FILE='piriformENDO96.LOTrast')
         WRITE (516,8789) time, L
          endif
        endif
       end do


       else IF (nodecell(thisno) .eq. 'semilunar   ') THEN
c semilunar  
c Determine which particular cells this node will be concerned with.
          firstcell = 1 
          lastcell =  num_semilunar 

        outrcd( 1) = time
        outrcd( 2) = v_semilunar   (1,firstcell+1)
        outrcd( 3) = v_semilunar   (numcomp_semilunar   ,firstcell+1)
c       outrcd( 3) = 0.01d0 * chi_semilunar   (48   ,firstcell+1)
        outrcd( 4) = v_semilunar   (48,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_semilunar(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_semilunar   
           z = z + gAMPA_semilunar   (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_semilunar   
           z = z + gNMDA_semilunar   (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_semilunar   
           z = z + gGABA_A_semilunar   (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
         z = 0.d0
          do i = 1, numcomp_semilunar   
           z = z + gGABA_B_semilunar   (i,firstcell+1)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 2 
        outrcd(10) = v_semilunar   (1, 226       )
        outrcd(11) = v_semilunar   (1, 439       )
        outrcd(12) = v_semilunar   (43,439        )
c       outrcd(12) = field_sup_tot   
c       outrcd(13) = field_deep_tot   
        outrcd(13) = 0.01d0 * chi_semilunar(48,firstcell+3)
        outrcd(14) = v_semilunar   (1,firstcell+4)
        outrcd(15) = v_semilunar   (numcomp_semilunar   ,firstcell+4)
          z = 0.d0
        do L = 1, num_semilunar
          if (ldistal_axon_semilunar(L).ge.0.d0) z = z + 1.d0
        end do
        outrcd(16) = z ! should be number of distal semilunar axons overshooting
        outrcd(17) = v_semilunar (1,firstcell + 450)
        outrcd(18) = v_semilunar (numcomp_semilunar,firstcell+450)
        outrcd(19) = v_semilunar (43,firstcell + 450)
c       outrcd(20) = v_semilunar (1,firstcell + 8)
        outrcd(20) = v_semilunar (1, 492         )
          if (place(thisno).eq.1) then
      OPEN(17,FILE='piriformENDO96.semilunar')
      WRITE (17,FMT='(20F11.3)') (OUTRCD(I),I=1,20)
          else
c           write(6,9091) 'semilunar', thisno, time, v_semilunar(1,firstcell),
c    &            v_semilunar(1,lastcell)
9091        format(a6,i4,3f10.4)
          endif



       do L = firstcell, lastcell     
        if (v_semilunar (1,L) .ge. -25.d0) then
          if (place(thisno).eq.1) then
         OPEN(416,FILE='piriformENDO96.semilunarrast')
         WRITE (416,8789) time, L
          endif
        endif
       end do

        outrcd( 1) = time
        outrcd( 2) = v_semilunar   (1,firstcell+3)
        outrcd( 3) = v_semilunar   (numcomp_semilunar   ,firstcell+3)
        outrcd( 4) = v_semilunar   (48,firstcell+3)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_semilunar(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_semilunar   
           z = z + gAMPA_semilunar   (i,firstcell+3)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_semilunar   
           z = z + gNMDA_semilunar   (i,firstcell+3)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_semilunar   
           z = z + gGABA_A_semilunar   (i,firstcell+3)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
         z = 0.d0
          do i = 1, numcomp_semilunar   
           z = z + gGABA_B_semilunar   (i,firstcell+3)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 2 
        outrcd(10) = v_semilunar   ( 2,firstcell+3)
        outrcd(11) = v_semilunar   ( 9,firstcell+3)
        outrcd(12) = v_semilunar   (43,firstcell+3)
        outrcd(13) = 0.01d0 * chi_semilunar(48,firstcell+3)
        outrcd(14) = v_semilunar   (31,firstcell+4)
        outrcd(15) = v_semilunar   ( 56              ,firstcell+3)
          z = 0.d0
        do L = 1, num_semilunar
          if (ldistal_axon_semilunar(L).ge.0.d0) z = z + 1.d0
        end do
        outrcd(16) = z ! should be number of distal semilunar axons overshooting
        outrcd(17) = v_semilunar (24,firstcell + 3)
        outrcd(18) = v_semilunar (25,firstcell + 3)
        outrcd(19) = v_semilunar (1, 492          )
        outrcd(20) = 1000.d0 * noisepe_semilunar 
          if (place(thisno).eq.1) then
c     OPEN(87,FILE='piriformENDO96.semilunarA')
c     WRITE (87,FMT='(20F10.4)') (OUTRCD(I),I=1,20)
c         else
          endif

       else IF (nodecell(thisno) .eq. 'multipolar') THEN
c multipolar  
c Determine which particular cells this node will be concerned with.
          firstcell = 1 
          lastcell =  num_multipolar 

        outrcd( 1) = time
        outrcd( 2) = v_multipolar   (1,firstcell+1)
        outrcd( 3) = v_multipolar (numcomp_multipolar  ,firstcell+1)
        outrcd( 4) = v_multipolar   (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_multipolar(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_multipolar   
           z = z + gAMPA_multipolar   (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_multipolar   
           z = z + gNMDA_multipolar   (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_multipolar   
           z = z + gGABA_A_multipolar   (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
         z = 0.d0
c         do i = 1, numcomp_multipolar   
c          z = z + gGABA_B_multipolar   (i,firstcell+1)
c         end do
c       outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 2 
         z = 0.d0
          do i = firstcell, lastcell
           if (V_multipolar (1,i).ge.0.d0) z = z+1.d0
          end do
        outrcd( 9) = z ! num overshooting somata 
        outrcd(10) = v_multipolar   (1,firstcell+2)
        outrcd(11) = v_multipolar   (1,firstcell+3)
        outrcd(12) = 0.d0 ! field_sup_tot   
        outrcd(13) = 0.d0 ! field_deep_tot    
        outrcd(14) = v_multipolar   (1,250        )
        outrcd(15) = v_multipolar   (1,lastcell-3 )
        outrcd(16) = v_multipolar   (1,lastcell-2 )
        outrcd(17) = v_multipolar   (1,lastcell-1 )
        outrcd(18) = v_multipolar   (1,lastcell   )
        outrcd(19) = curr_multipolar(1,1) 
        outrcd(20) = 1.d3 * gapcon_multipolar 
          if (place(thisno).eq.1) then
      OPEN(18,FILE='piriformENDO96.multipolar')
      WRITE (18,FMT='(20F10.4)') (OUTRCD(I),I=1,20)
          end if

       do L = firstcell, lastcell     
        if (v_multipolar (1,L) .ge. -10.d0) then
         OPEN(621,FILE='piriformENDO96.multipolarrast')
         WRITE (621,8789) time, L
        endif
       end do

       else IF (nodecell(thisno) .eq. 'L3pyr       ') THEN
c L3pyr
c Determine which particular cells this node will be concerned with.
         if (mod(O,500).eq.0) then
          write(6,3835) time, v_L3pyr    (1,5)
3835      format(' L3pyr     ',f10.4,2x,f10.3)
         endif ! just to make sure job is running

          firstcell = 1 
          lastcell =  num_L3pyr 

        outrcd( 1) = time
        outrcd( 2) = v_L3pyr(1,firstcell+1)
        outrcd( 3) = v_L3pyr(numcomp_L3pyr,firstcell+1)
        outrcd( 4) = v_L3pyr(43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_L3pyr(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_L3pyr
           z = z + gAMPA_L3pyr(i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_L3pyr
           z = z + gNMDA_L3pyr(i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_L3pyr
           z = z + gGABA_A_L3pyr(i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
         z = 0.d0
          do i = 1, numcomp_L3pyr
           z = z + gGABA_B_L3pyr(i,firstcell+1)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 2 
        outrcd(10) = v_L3pyr(1, 96        )
        outrcd(11) = v_L3pyr(1,firstcell+3)
        outrcd(12) = 0.d0 ! field_sup_tot      
        outrcd(13) = 0.d0 ! field_deep_tot       
        outrcd(14) = v_L3pyr(1,  22           )
        outrcd(15) = v_L3pyr(43, 22           )
        outrcd(16) = v_L3pyr(1, 388         )
          if (place(thisno).eq.1) then
      OPEN(19,FILE='piriformENDO96.L3pyr')
      WRITE (19,FMT='(16F10.4)') (OUTRCD(I),I=1,16)
c         else
c      write(6,9092) 'L3pyr',thisno,time,v_L3pyr(1,firstcell),
c    &            v_L3pyr(1,lastcell)
c9092        format(a9,i4,3f10.4)
          endif

       do L = firstcell, lastcell     
        if (v_L3pyr (1,L) .ge. -25.d0) then
          if (place(thisno).eq.1) then
         OPEN(415,FILE='piriformENDO96.L3pyrrast')
         WRITE (415,8789) time, L
          endif
        end if
       end do

       else IF (nodecell(thisno) .eq. 'deepintern  ') THEN
c deepbask 
c Determine which particular cells this node will be concerned with.
          firstcell = 1 
          lastcell =  num_deepbask 

        outrcd( 1) = time
        outrcd( 2) = v_deepbask (1,firstcell+1)
        outrcd( 3) = v_deepbask (numcomp_deepbask ,firstcell+1)
        outrcd( 4) = v_deepbask (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_deepbask(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_deepbask 
           z = z + gAMPA_deepbask (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepbask 
           z = z + gNMDA_deepbask (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepbask 
           z = z + gGABA_A_deepbask (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_deepbask (1,firstcell+2)
        outrcd(10) = v_deepbask (1,firstcell+3)
           if (place(thisno).eq.1) then
      OPEN(20,FILE='piriformENDO96.deepbask')
      WRITE (20,FMT='(10F10.4)') (OUTRCD(I),I=1,10)
           endif

       do L = firstcell, lastcell     
        if (v_deepbask  (1,L) .ge.   0.d0) then
         OPEN(516,FILE='piriformENDO96.deepbaskrast')
         WRITE (516,8789) time, L
        endif
       end do

c deepng 
c Determine which particular cells this node will be concerned with.
          firstcell = 1 
          lastcell =  num_deepng 

        outrcd( 1) = time
        outrcd( 2) = v_deepng (1,firstcell+1)
        outrcd( 3) = v_deepng (numcomp_deepng ,firstcell+1)
        outrcd( 4) = v_deepng (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_deepng(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_deepng 
           z = z + gAMPA_deepng (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepng 
           z = z + gNMDA_deepng (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepng 
           z = z + gGABA_A_deepng (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_deepng (1,firstcell+2)
        outrcd(10) = v_deepng (1,firstcell+3)
      OPEN(34,FILE='piriformENDO96.deepng')
      WRITE (34,FMT='(10F10.4)') (OUTRCD(I),I=1,10)


c deepLTS
c Determine which particular cells this node will be concerned with.
          firstcell = 1 
          lastcell =  num_deepLTS 

        outrcd( 1) = time
        outrcd( 2) = v_deepLTS (1,firstcell+1)
        outrcd( 3) = v_deepLTS (numcomp_deepLTS ,firstcell+1)
        outrcd( 4) = v_deepLTS (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_deepLTS(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_deepLTS 
           z = z + gAMPA_deepLTS (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepLTS 
           z = z + gNMDA_deepLTS (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_deepLTS 
           z = z + gGABA_A_deepLTS (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
        outrcd( 9) = v_deepLTS (1,firstcell+2)
        outrcd(10) = v_deepLTS (1,firstcell+3)
           if (place(thisno).eq.1) then
      OPEN(21,FILE='piriformENDO96.deepLTS')
      WRITE (21,FMT='(10F10.4)') (OUTRCD(I),I=1,10)
           endif

       do L = firstcell, lastcell     
        if (v_deepLTS   (1,L) .ge.   0.d0) then
         OPEN(616,FILE='piriformENDO96.deepLTSrast')
         WRITE (616,8789) time, L
        endif
       end do

       else IF (nodecell(thisno) .eq. 'placeholder5') THEN
c placeholder5      
c Determine which particular cells this node will be concerned with.
          firstcell = 1 
          lastcell =  num_placeholder5 

        outrcd( 1) = time
        outrcd( 2) = v_placeholder5      (1,firstcell+1)
        outrcd( 3) = v_placeholder5  (numcomp_placeholder5 ,firstcell+1)
        outrcd( 4) = v_placeholder5      (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_placeholder5(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_placeholder5      
           z = z + gAMPA_placeholder5      (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_placeholder5      
           z = z + gNMDA_placeholder5      (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_placeholder5      
           z = z + gGABA_A_placeholder5      (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
         z = 0.d0
          do i = 1, numcomp_placeholder5      
           z = z + gGABA_B_placeholder5      (i,firstcell+1)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 2 
        outrcd(10) = v_placeholder5      (1,firstcell+2)
        outrcd(11) = v_placeholder5      (1,firstcell+3)

          z = 0.d0
          do i = firstcell, lastcell
           if(v_placeholder5 (numcomp_placeholder5,i).gt.0.d0) z = z + 1.d0
          end do
        outrcd(12) = z   

        outrcd(13) = v_placeholder5 (1,33)

c     OPEN(23,FILE='piriformENDO96.placeholder5')
c     WRITE (23,FMT='(13F10.4)') (OUTRCD(I),I=1,13)

       else IF (nodecell(thisno) .eq. 'placeholder6') THEN
c placeholder6       
c Determine which particular cells this node will be concerned with.
          firstcell = 1 
          lastcell =  num_placeholder6 

        outrcd( 1) = time
        outrcd( 2) = v_placeholder6      (1,firstcell+1)
        outrcd( 3) = v_placeholder6  (numcomp_placeholder6 ,firstcell+1)
        outrcd( 4) = v_placeholder6      (43,firstcell+1)
         z = 0.d0
          do i = firstcell, lastcell
           z = z - v_placeholder6(1,i)
          end do
        outrcd( 5) = z / dble(lastcell-firstcell+1) !  -av. cell somata 
         z = 0.d0
          do i = 1, numcomp_placeholder6       
           z = z + gAMPA_placeholder6      (i,firstcell+1)
          end do
        outrcd( 6) = z * 1000.d0 ! total AMPA cell 2 
         z = 0.d0
          do i = 1, numcomp_placeholder6         
           z = z + gNMDA_placeholder6      (i,firstcell+1)
          end do
        outrcd( 7) = z * 1000.d0 ! total NMDA cell 2 
         z = 0.d0
          do i = 1, numcomp_placeholder6        
           z = z + gGABA_A_placeholder6      (i,firstcell+1)
          end do
        outrcd( 8) = z * 1000.d0 ! total GABA-A, cell 2 
         z = 0.d0
          do i = 1, numcomp_placeholder6        
           z = z + gGABA_B_placeholder6      (i,firstcell+1)
          end do
        outrcd( 9) = z * 1000.d0 ! total GABA-B, cell 2 
        outrcd(10) = v_placeholder6      (1,firstcell+2)
        outrcd(11) = v_placeholder6      (1,firstcell+3)

          z = 0.d0
          do i = firstcell, lastcell
           if(v_placeholder6 (numcomp_placeholder6 ,i).gt.0.d0) z= z + 1.d0
          end do
        outrcd(12) = z   
c     OPEN(24,FILE='piriformENDO96.placeholder6')
c     WRITE (24,FMT='(12F10.4)') (OUTRCD(I),I=1,12)
       endif ! checking thisno

       endif ! mod(O, 50) = 0

        goto 1000
c END guts of main program

2000    CONTINUE
        time2 = gettime()
         if (thisno.eq.0) then
        write(6,3434) time2 - time1
         endif
3434    format(' elapsed time = ',f8.0,' secs')

        call mpi_finalize (info)
             END
c end main routine


c 22 Aug 2019, start with suppyrRS integration subroutine from
c son_of_groucho, and use for L2pyr in piriform simulations.
c Need to change field variables and depth definitions,
c  and perhaps alter compartment dimensions.
c 11 Sept 2006, start with /interact/integrate_suppyrRSXP.f & add GABA-B
! 7 Nov. 2005: modify integrate_suppyrRSX.f to allow for Colbert-Pan axon.
!29 July 2005: modify groucho/integrate_suppyrRS.f, for a separate
! call for initialization, and to integrate only selected cells.
! Integration routine for suppyrRS cells
! Routine adapted from scortn in supergj.f
c      SUBROUTINE INTEGRATE_suppyrRSXPB (O, time, numcell,     
       SUBROUTINE INTEGRATE_L2pyr       (O, time, numcell,     
     &    V, curr, initialize, firstcell, lastcell,
     & gAMPA, gNMDA, gGABA_A, gGABA_B,
     & Mg, 
     & gapcon  ,totaxgj   ,gjtable, dt,
     &  chi,mnaf,mnap,
     &  hnaf,mkdr,mka,
     &  hka,mk2,hk2,
     &  mkm,mkc,mkahp,
     &  mcat,hcat,mcal,
     &  mar,field_sup,field_deep,rel_axonshift)

       SAVE

       INTEGER, PARAMETER:: numcomp = 74
! numcomp here must be compatible with numcomp_suppyrRS in calling prog.
       INTEGER  numcell, num_other
       INTEGER initialize, firstcell, lastcell
       INTEGER J1, I, J, K, K1, K2, K3, L, L1, O
       REAL*8 c(numcomp), curr(numcomp,numcell)
       REAL*8  Z, Z1, Z2, Z3, Z4, DT, time
       integer totaxgj, gjtable(totaxgj,4)
       real*8 gapcon, gAMPA(numcomp,numcell),
     &        gNMDA(numcomp,numcell), gGABA_A(numcomp,numcell),
     &        gGABA_B(numcomp,numcell)
       real*8 Mg, V(numcomp,numcell), rel_axonshift

c CINV is 1/C, i.e. inverse capacitance
       real*8 chi(numcomp,numcell),
     & mnaf(numcomp,numcell),mnap(numcomp,numcell),
     x hnaf(numcomp,numcell), mkdr(numcomp,numcell),
     x mka(numcomp,numcell),hka(numcomp,numcell),
     x mk2(numcomp,numcell), cinv(numcomp),
     x hk2(numcomp,numcell),mkm(numcomp,numcell),
     x mkc(numcomp,numcell),mkahp(numcomp,numcell),
     x mcat(numcomp,numcell),hcat(numcomp,numcell),
     x mcal(numcomp,numcell), betchi(numcomp),
     x mar(numcomp,numcell),jacob(numcomp,numcomp),
     x gam(0: numcomp,0: numcomp),gL(numcomp),gnaf(numcomp),
     x gnap(numcomp),gkdr(numcomp),gka(numcomp),
     x gk2(numcomp),gkm(numcomp),
     x gkc(numcomp),gkahp(numcomp),
     x gcat(numcomp),gcaL(numcomp),gar(numcomp),
     x cafor(numcomp)
       real*8
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
       real*8 vL(numcomp),vk(numcomp),vna,var,vca,vgaba_a
       real*8 depth(12), membcurr(12), field_sup, field_deep
       integer level(numcomp)

        INTEGER NEIGH(numcomp,10), NNUM(numcomp)
        INTEGER igap1, igap2
c the f's are the functions giving 1st derivatives for evolution of
c the differential equations for the voltages (v), calcium (chi), and
c other state variables.
       real*8 fv(numcomp), fchi(numcomp),
     x fmnaf(numcomp),fhnaf(numcomp),fmkdr(numcomp),
     x fmka(numcomp),fhka(numcomp),fmk2(numcomp),
     x fhk2(numcomp),fmnap(numcomp),
     x fmkm(numcomp),fmkc(numcomp),fmkahp(numcomp),
     x fmcat(numcomp),fhcat(numcomp),fmcal(numcomp),
     x fmar(numcomp)

c below are for calculating the partial derivatives
       real*8 dfv_dv(numcomp,numcomp), dfv_dchi(numcomp),
     x  dfv_dmnaf(numcomp),  dfv_dmnap(numcomp),
     x  dfv_dhnaf(numcomp),dfv_dmkdr(numcomp),
     x  dfv_dmka(numcomp),dfv_dhka(numcomp),
     x  dfv_dmk2(numcomp),dfv_dhk2(numcomp),
     x  dfv_dmkm(numcomp),dfv_dmkc(numcomp),
     x  dfv_dmkahp(numcomp),dfv_dmcat(numcomp),
     x  dfv_dhcat(numcomp),dfv_dmcal(numcomp),
     x  dfv_dmar(numcomp)

        real*8 dfchi_dv(numcomp), dfchi_dchi(numcomp),
     x dfmnaf_dmnaf(numcomp), dfmnaf_dv(numcomp),
     x dfhnaf_dhnaf(numcomp),
     x dfmnap_dmnap(numcomp), dfmnap_dv(numcomp),
     x dfhnaf_dv(numcomp),dfmkdr_dmkdr(numcomp),
     x dfmkdr_dv(numcomp),
     x dfmka_dmka(numcomp),dfmka_dv(numcomp),
     x dfhka_dhka(numcomp),dfhka_dv(numcomp),
     x dfmk2_dmk2(numcomp),dfmk2_dv(numcomp),
     x dfhk2_dhk2(numcomp),dfhk2_dv(numcomp),
     x dfmkm_dmkm(numcomp),dfmkm_dv(numcomp),
     x dfmkc_dmkc(numcomp),dfmkc_dv(numcomp),
     x dfmcat_dmcat(numcomp),dfmcat_dv(numcomp),dfhcat_dhcat(numcomp),
     x dfhcat_dv(numcomp),dfmcal_dmcal(numcomp),dfmcal_dv(numcomp),
     x dfmar_dmar(numcomp),dfmar_dv(numcomp),dfmkahp_dchi(numcomp),
     x dfmkahp_dmkahp(numcomp), dt2

       REAL*8 OPEN(numcomp),gamma(numcomp),gamma_prime(numcomp)
c gamma is function of chi used in calculating KC conductance
       REAL*8 alpham_ahp(numcomp), alpham_ahp_prime(numcomp)
       REAL*8 gna_tot(numcomp),gk_tot(numcomp),gca_tot(numcomp)
       REAL*8 gca_high(numcomp), gar_tot(numcomp)
c this will be gCa conductance corresponding to high-thresh channels

       real*8 persistentNa_shift, fastNa_shift_SD,
     x   fastNa_shift_axon

       REAL*8 A, BB1, BB2  ! params. for FNMDA.f


c          if (O.eq.1) then
           if (initialize.eq.0) then
c do initialization

c Program fnmda assumes A, BB1, BB2 defined in calling program
c as follows:
         A = DEXP(-2.847d0)
         BB1 = DEXP(-.693d0)
         BB2 = DEXP(-3.101d0)

c       goto 4000
       CALL   SCORT_SETUP_L2pyr   
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)

        CALL SCORTMAJ_L2pyr   
     X             (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM,depth,level)

          do i = 1, numcomp
             cinv(i) = 1.d0 / c(i)
          end do
4000      CONTINUE

           do i = 1, numcomp
          vL(i) = -70.d0
          vK(i) = -95.d0
           end do

        VNA = 50.d0
        VCA = 125.d0
        VAR = -43.d0
        VAR = -35.d0
c -43 mV from Huguenard & McCormick
        VGABA_A = -81.d0
c       write(6,901) VNa, VCa, VK(1), O
901     format('VNa =',f6.2,' VCa =',f6.2,' VK =',f6.2,
     &   ' O = ',i3)

c ? initialize membrane state variables?
         do L = 1, numcell  
         do i = 1, numcomp
        v(i,L) = VL(i)
	chi(i,L) = 0.d0
	mnaf(i,L) = 0.d0
	mkdr(i,L) = 0.d0
	mk2(i,L) = 0.d0
	mkm(i,L) = 0.d0
	mkc(i,L) = 0.d0
	mkahp(i,L) = 0.d0
	mcat(i,L) = 0.d0
	mcal(i,L) = 0.d0
         end do
         end do

          do L = 1, numcell
        k1 = idnint (4.d0 * (v(1,L) + 120.d0))

            do i = 1, numcomp
      hnaf(i,L) = alphah_naf(k1)/(alphah_naf(k1)
     &       +betah_naf(k1))
      hka(i,L) = alphah_ka(k1)/(alphah_ka(k1)
     &                               +betah_ka(k1))
      hk2(i,L) = alphah_k2(k1)/(alphah_k2(k1)
     &                                +betah_k2(k1))
      hcat(i,L)=alphah_cat(k1)/(alphah_cat(k1)
     &                                +betah_cat(k1))
c     mar=alpham_ar(k1)/(alpham_ar(k1)+betam_ar(k1))
      mar(i,L) = .25d0
             end do
           end do


             do i = 1, numcomp
	    open(i) = 0.d0
            gkm(i) = 2.d0 * gkm(i)
             end do

         do i = 1, 68
c          gnaf(i) = 0.8d0 * 1.25d0 * gnaf(i) ! factor of 0.8 added 19 Nov. 2005
c          gnaf(i) = 0.9d0 * 1.25d0 * gnaf(i) ! Back to 0.9, 29 Nov. 2005
           gnaf(i) = 0.6d0 * 1.25d0 * gnaf(i) ! 
! NOTE THAT THERE IS QUESTION OF HOW TO COMPARE BEHAVIOR OF PYRAMID IN NETWORK WITH
! SIMULATIONS OF SINGLE CELL.  IN FORMER CASE, THERE IS LARGE AXONAL SHUNT THROUGH
! gj(s), NOT PRESENT IN SINGLE CELL MODEL.  THEREFORE, HIGHER AXONAL gNa MIGHT BE
! NECESSARY FOR SPIKE PROPAGATION.
c          gnaf(i) = 0.9d0 * 1.25d0 * gnaf(i) ! factor of 0.9 added 20 Nov. 2005
           gkdr(i) = 1.25d0 * gkdr(i)
         end do
 
c Perhaps reduce fast gNa on IS
          gnaf(69) = 1.00d0 * gnaf(69)
c         gnaf(69) = 0.25d0 * gnaf(69)
          gnaf(70) = 1.00d0 * gnaf(70)
c         gnaf(70) = 0.25d0 * gnaf(70)

c Perhaps reduce coupling between soma and IS
c         gam(1,69) = 0.15d0 * gam(1,69)
c         gam(69,1) = 0.15d0 * gam(69,1)

               z1 = 0.0d0
c              z2 = 1.2d0 ! value 1.2 tried Feb. 21, 2013
c              z2 = 1.5d0 ! value 1.2 tried Feb. 21, 2013
               z2 = 2.0d0 ! value 1.2 tried Feb. 21, 2013
               z3 = 1.0d0
c              z3 = 0.4d0
c              z3 = 0.0d0 ! Note reduction from 0.4, to prevent
c slow hyperpolarization that seems to mess up gamma.
               z4 = 0.3d0
c RS cell
             do i = 1, numcomp
              gnap(i) = z1 * gnap(i)
              gkc (i) = z2 * gkc (i)
              gkahp(i) = z3 * gkahp(i)
              gkm (i) = z4 * gkm (i)
             end do

              goto 6000

          endif
c End initialization

          do i = 1, 12
           membcurr(i) = 0.d0
          end do

c                  goto 2001


c             do L = 1, numcell
              do L = firstcell, lastcell

	  do i = 1, numcomp
	  do j = 1, nnum(i)
	   if (neigh(i,j).gt.numcomp) then
          write(6,433) i, j, L
433       format(' ls ',3x,3i5)
           endif
	end do
	end do

       DO I = 1, numcomp
          FV(I) = -GL(I) * (V(I,L) - VL(i)) * cinv(i)
          DO J = 1, NNUM(I)
             K = NEIGH(I,J)
302     FV(I) = FV(I) + GAM(I,K) * (V(K,L) - V(I,L)) * cinv(i)
           END DO
       END DO
301    CONTINUE


       CALL FNMDA (V, OPEN, numcell, numcomp, MG, L,
     &                 A, BB1, BB2)

      DO I = 1, numcomp
       FV(I) = FV(I) + ( CURR(I,L)
     X   - (gampa(I,L) + open(i) * gnmda(I,L))*V(I,L)
     X   - ggaba_a(I,L)*(V(I,L)-Vgaba_a) 
     X   - ggaba_b(I,L)*(V(I,L)-VK(i)  ) ) * cinv(i)
c above assumes equil. potential for AMPA & NMDA = 0 mV
      END DO
421      continue

       do m = 1, totaxgj
        if (gjtable(m,1).eq.L) then
         L1 = gjtable(m,3)
         igap1 = gjtable(m,2)
         igap2 = gjtable(m,4)
 	fv(igap1) = fv(igap1) + gapcon *
     &   (v(igap2,L1) - v(igap1,L)) * cinv(igap1)
        else if (gjtable(m,3).eq.L) then
         L1 = gjtable(m,1)
         igap1 = gjtable(m,4)
         igap2 = gjtable(m,2)
 	fv(igap1) = fv(igap1) + gapcon *
     &   (v(igap2,L1) - v(igap1,L)) * cinv(igap1)
        endif
       end do ! do m


       do i = 1, numcomp
        gamma(i) = dmin1 (1.d0, .004d0 * chi(i,L))
        if (chi(i,L).le.250.d0) then
          gamma_prime(i) = .004d0
        else
          gamma_prime(i) = 0.d0
        endif
c         endif
       end do

      DO I = 1, numcomp
       gna_tot(i) = gnaf(i) * (mnaf(i,L)**3) * hnaf(i,L) +
     x     gnap(i) * mnap(i,L)
       gk_tot(i) = gkdr(i) * (mkdr(i,L)**4) +
     x             gka(i)  * (mka(i,L)**4) * hka(i,L) +
     x             gk2(i)  * mk2(i,L) * hk2(i,L) +
     x             gkm(i)  * mkm(i,L) +
     x             gkc(i)  * mkc(i,L) * gamma(i) +
     x             gkahp(i)* mkahp(i,L)
       gca_tot(i) = gcat(i) * (mcat(i,L)**2) * hcat(i,L) +
     x              gcaL(i) * (mcaL(i,L)**2)
       gca_high(i) =
     x              gcaL(i) * (mcaL(i,L)**2)
       gar_tot(i) = gar(i) * mar(i,L)


       FV(I) = FV(I) - ( gna_tot(i) * (v(i,L) - vna)
     X  + gk_tot(i) * (v(i,L) - vK(i))
     X  + gca_tot(i) * (v(i,L) - vCa)
     X  + gar_tot(i) * (v(i,L) - var) ) * cinv(i)
c        endif
       END DO
88           continue

         do i = 1, numcomp
         do j = 1, numcomp
          if (i.ne.j) then
            dfv_dv(i,j) = jacob(i,j)
          else
            dfv_dv(i,j) = jacob(i,i) - cinv(i) *
     X  (gna_tot(i) + gk_tot(i) + gca_tot(i) + gar_tot(i)
     X   + ggaba_a(i,L) + ggaba_b(i,L) + gampa(i,L)
     X   + open(i) * gnmda(I,L) )
          endif
         end do
         end do

           do i = 1, numcomp
        dfv_dchi(i)  = - cinv(i) * gkc(i) * mkc(i,L) *
     x                     gamma_prime(i) * (v(i,L)-vK(i))
        dfv_dmnaf(i) = -3.d0 * cinv(i) * (mnaf(i,L)**2) *
     X    (gnaf(i) * hnaf(i,L)          ) * (v(i,L) - vna)
        dfv_dmnap(i) = - cinv(i) *
     X    (               gnap(i)) * (v(i,L) - vna)
        dfv_dhnaf(i) = - cinv(i) * gnaf(i) * (mnaf(i,L)**3) *
     X                    (v(i,L) - vna)
        dfv_dmkdr(i) = -4.d0 * cinv(i) * gkdr(i) * (mkdr(i,L)**3)
     X                   * (v(i,L) - vK(i))
        dfv_dmka(i)  = -4.d0 * cinv(i) * gka(i) * (mka(i,L)**3) *
     X                   hka(i,L) * (v(i,L) - vK(i))
        dfv_dhka(i)  = - cinv(i) * gka(i) * (mka(i,L)**4) *
     X                    (v(i,L) - vK(i))
      dfv_dmk2(i) = - cinv(i) * gk2(i) * hk2(i,L) * (v(i,L)-vK(i))
      dfv_dhk2(i) = - cinv(i) * gk2(i) * mk2(i,L) * (v(i,L)-vK(i))
      dfv_dmkm(i) = - cinv(i) * gkm(i) * (v(i,L) - vK(i))
      dfv_dmkc(i) = - cinv(i)*gkc(i) * gamma(i) * (v(i,L)-vK(i))
        dfv_dmkahp(i)= - cinv(i) * gkahp(i) * (v(i,L) - vK(i))
        dfv_dmcat(i)  = -2.d0 * cinv(i) * gcat(i) * mcat(i,L) *
     X                    hcat(i,L) * (v(i,L) - vCa)
        dfv_dhcat(i) = - cinv(i) * gcat(i) * (mcat(i,L)**2) *
     X                  (v(i,L) - vCa)
        dfv_dmcal(i) = -2.d0 * cinv(i) * gcal(i) * mcal(i,L) *
     X                      (v(i,L) - vCa)
        dfv_dmar(i) = - cinv(i) * gar(i) * (v(i,L) - var)
            end do

         do i = 1, numcomp
          fchi(i) = - cafor(i) * gca_high(i) * (v(i,L) - vca)
     x       - betchi(i) * chi(i,L)
          dfchi_dv(i) = - cafor(i) * gca_high(i)
          dfchi_dchi(i) = - betchi(i)
         end do

       do i = 1, numcomp
c Note possible increase in rate at which AHP current develops
c       alpham_ahp(i) = dmin1(0.2d-4 * chi(i,L),0.01d0)
        alpham_ahp(i) = dmin1(1.0d-4 * chi(i,L),0.01d0)
        if (chi(i,L).le.500.d0) then
c         alpham_ahp_prime(i) = 0.2d-4
          alpham_ahp_prime(i) = 1.0d-4
        else
          alpham_ahp_prime(i) = 0.d0
        endif
       end do

       do i = 1, numcomp
        fmkahp(i) = alpham_ahp(i) * (1.d0 - mkahp(i,L))
c    x                  -.001d0 * mkahp(i,L)
     x                  -.010d0 * mkahp(i,L)
c       dfmkahp_dmkahp(i) = - alpham_ahp(i) - .001d0
        dfmkahp_dmkahp(i) = - alpham_ahp(i) - .010d0
        dfmkahp_dchi(i) = alpham_ahp_prime(i) *
     x                     (1.d0 - mkahp(i,L))
       end do

          do i = 1, numcomp

       K1 = IDNINT ( 4.d0 * (V(I,L) + 120.d0) )
       IF (K1.GT.640) K1 = 640
       IF (K1.LT.  0) K1 =   0

c      persistentNa_shift =  0.d0
c      persistentNa_shift =  8.d0
       persistentNa_shift = 10.d0
       K2 = IDNINT ( 4.d0 * (V(I,L)+persistentNa_shift+ 120.d0) )
       IF (K2.GT.640) K2 = 640
       IF (K2.LT.  0) K2 =   0

c            fastNa_shift = -2.0d0
c            fastNa_shift = -2.5d0
             fastNa_shift_SD = -3.5d0
             fastNa_shift_axon = fastNa_shift_SD + rel_axonshift 
       K0 = IDNINT ( 4.d0 * (V(I,L)+  fastNa_shift_SD+ 120.d0) )
       IF (K0.GT.640) K0 = 640
       IF (K0.LT.  0) K0 =   0
       K3 = IDNINT ( 4.d0 * (V(I,L)+  fastNa_shift_axon+ 120.d0) )
       IF (K3.GT.640) K3 = 640
       IF (K3.LT.  0) K3 =   0

         if (i.le.68) then   ! FOR SD
        fmnaf(i) = alpham_naf(k0) * (1.d0 - mnaf(i,L)) -
     X              betam_naf(k0) * mnaf(i,L)
        fhnaf(i) = alphah_naf(k0) * (1.d0 - hnaf(i,L)) -
     X              betah_naf(k0) * hnaf(i,L)
         else  ! for axon
        fmnaf(i) = alpham_naf(k3) * (1.d0 - mnaf(i,L)) -
     X              betam_naf(k3) * mnaf(i,L)
        fhnaf(i) = alphah_naf(k3) * (1.d0 - hnaf(i,L)) -
     X              betah_naf(k3) * hnaf(i,L)
         endif
        fmnap(i) = alpham_naf(k2) * (1.d0 - mnap(i,L)) -
     X              betam_naf(k2) * mnap(i,L)
        fmkdr(i) = alpham_kdr(k1) * (1.d0 - mkdr(i,L)) -
     X              betam_kdr(k1) * mkdr(i,L)
        fmka(i)  = alpham_ka (k1) * (1.d0 - mka(i,L)) -
     X              betam_ka (k1) * mka(i,L)
        fhka(i)  = alphah_ka (k1) * (1.d0 - hka(i,L)) -
     X              betah_ka (k1) * hka(i,L)
        fmk2(i)  = alpham_k2 (k1) * (1.d0 - mk2(i,L)) -
     X              betam_k2 (k1) * mk2(i,L)
        fhk2(i)  = alphah_k2 (k1) * (1.d0 - hk2(i,L)) -
     X              betah_k2 (k1) * hk2(i,L)
        fmkm(i)  = alpham_km (k1) * (1.d0 - mkm(i,L)) -
     X              betam_km (k1) * mkm(i,L)
        fmkc(i)  = alpham_kc (k1) * (1.d0 - mkc(i,L)) -
     X              betam_kc (k1) * mkc(i,L)
        fmcat(i) = alpham_cat(k1) * (1.d0 - mcat(i,L)) -
     X              betam_cat(k1) * mcat(i,L)
        fhcat(i) = alphah_cat(k1) * (1.d0 - hcat(i,L)) -
     X              betah_cat(k1) * hcat(i,L)
        fmcaL(i) = alpham_caL(k1) * (1.d0 - mcaL(i,L)) -
     X              betam_caL(k1) * mcaL(i,L)
        fmar(i)  = alpham_ar (k1) * (1.d0 - mar(i,L)) -
     X              betam_ar (k1) * mar(i,L)

       dfmnaf_dv(i) = dalpham_naf(k0) * (1.d0 - mnaf(i,L)) -
     X                  dbetam_naf(k0) * mnaf(i,L)
       dfmnap_dv(i) = dalpham_naf(k2) * (1.d0 - mnap(i,L)) -
     X                  dbetam_naf(k2) * mnap(i,L)
       dfhnaf_dv(i) = dalphah_naf(k1) * (1.d0 - hnaf(i,L)) -
     X                  dbetah_naf(k1) * hnaf(i,L)
       dfmkdr_dv(i) = dalpham_kdr(k1) * (1.d0 - mkdr(i,L)) -
     X                  dbetam_kdr(k1) * mkdr(i,L)
       dfmka_dv(i)  = dalpham_ka(k1) * (1.d0 - mka(i,L)) -
     X                  dbetam_ka(k1) * mka(i,L)
       dfhka_dv(i)  = dalphah_ka(k1) * (1.d0 - hka(i,L)) -
     X                  dbetah_ka(k1) * hka(i,L)
       dfmk2_dv(i)  = dalpham_k2(k1) * (1.d0 - mk2(i,L)) -
     X                  dbetam_k2(k1) * mk2(i,L)
       dfhk2_dv(i)  = dalphah_k2(k1) * (1.d0 - hk2(i,L)) -
     X                  dbetah_k2(k1) * hk2(i,L)
       dfmkm_dv(i)  = dalpham_km(k1) * (1.d0 - mkm(i,L)) -
     X                  dbetam_km(k1) * mkm(i,L)
       dfmkc_dv(i)  = dalpham_kc(k1) * (1.d0 - mkc(i,L)) -
     X                  dbetam_kc(k1) * mkc(i,L)
       dfmcat_dv(i) = dalpham_cat(k1) * (1.d0 - mcat(i,L)) -
     X                  dbetam_cat(k1) * mcat(i,L)
       dfhcat_dv(i) = dalphah_cat(k1) * (1.d0 - hcat(i,L)) -
     X                  dbetah_cat(k1) * hcat(i,L)
       dfmcaL_dv(i) = dalpham_caL(k1) * (1.d0 - mcaL(i,L)) -
     X                  dbetam_caL(k1) * mcaL(i,L)
       dfmar_dv(i)  = dalpham_ar(k1) * (1.d0 - mar(i,L)) -
     X                  dbetam_ar(k1) * mar(i,L)

       dfmnaf_dmnaf(i) =  - alpham_naf(k0) - betam_naf(k0)
       dfmnap_dmnap(i) =  - alpham_naf(k2) - betam_naf(k2)
       dfhnaf_dhnaf(i) =  - alphah_naf(k1) - betah_naf(k1)
       dfmkdr_dmkdr(i) =  - alpham_kdr(k1) - betam_kdr(k1)
       dfmka_dmka(i)  =   - alpham_ka (k1) - betam_ka (k1)
       dfhka_dhka(i)  =   - alphah_ka (k1) - betah_ka (k1)
       dfmk2_dmk2(i)  =   - alpham_k2 (k1) - betam_k2 (k1)
       dfhk2_dhk2(i)  =   - alphah_k2 (k1) - betah_k2 (k1)
       dfmkm_dmkm(i)  =   - alpham_km (k1) - betam_km (k1)
       dfmkc_dmkc(i)  =   - alpham_kc (k1) - betam_kc (k1)
       dfmcat_dmcat(i) =  - alpham_cat(k1) - betam_cat(k1)
       dfhcat_dhcat(i) =  - alphah_cat(k1) - betah_cat(k1)
       dfmcaL_dmcaL(i) =  - alpham_caL(k1) - betam_caL(k1)
       dfmar_dmar(i)  =   - alpham_ar (k1) - betam_ar (k1)

          end do

       dt2 = 0.5d0 * dt * dt

        do i = 1, numcomp
          v(i,L) = v(i,L) + dt * fv(i)
           do j = 1, numcomp
        v(i,L) = v(i,L) + dt2 * dfv_dv(i,j) * fv(j)
           end do
        v(i,L) = v(i,L) + dt2 * ( dfv_dchi(i) * fchi(i)
     X          + dfv_dmnaf(i) * fmnaf(i)
     X          + dfv_dmnap(i) * fmnap(i)
     X          + dfv_dhnaf(i) * fhnaf(i)
     X          + dfv_dmkdr(i) * fmkdr(i)
     X          + dfv_dmka(i)  * fmka(i)
     X          + dfv_dhka(i)  * fhka(i)
     X          + dfv_dmk2(i)  * fmk2(i)
     X          + dfv_dhk2(i)  * fhk2(i)
     X          + dfv_dmkm(i)  * fmkm(i)
     X          + dfv_dmkc(i)  * fmkc(i)
     X          + dfv_dmkahp(i)* fmkahp(i)
     X          + dfv_dmcat(i)  * fmcat(i)
     X          + dfv_dhcat(i) * fhcat(i)
     X          + dfv_dmcaL(i) * fmcaL(i)
     X          + dfv_dmar(i)  * fmar(i) )

        chi(i,L) = chi(i,L) + dt * fchi(i) + dt2 *
     X   (dfchi_dchi(i) * fchi(i) + dfchi_dv(i) * fv(i))
        mnaf(i,L) = mnaf(i,L) + dt * fmnaf(i) + dt2 *
     X   (dfmnaf_dmnaf(i) * fmnaf(i) + dfmnaf_dv(i)*fv(i))
        mnap(i,L) = mnap(i,L) + dt * fmnap(i) + dt2 *
     X   (dfmnap_dmnap(i) * fmnap(i) + dfmnap_dv(i)*fv(i))
        hnaf(i,L) = hnaf(i,L) + dt * fhnaf(i) + dt2 *
     X   (dfhnaf_dhnaf(i) * fhnaf(i) + dfhnaf_dv(i)*fv(i))
        mkdr(i,L) = mkdr(i,L) + dt * fmkdr(i) + dt2 *
     X   (dfmkdr_dmkdr(i) * fmkdr(i) + dfmkdr_dv(i)*fv(i))
        mka(i,L) =  mka(i,L) + dt * fmka(i) + dt2 *
     X   (dfmka_dmka(i) * fmka(i) + dfmka_dv(i) * fv(i))
        hka(i,L) =  hka(i,L) + dt * fhka(i) + dt2 *
     X   (dfhka_dhka(i) * fhka(i) + dfhka_dv(i) * fv(i))
        mk2(i,L) =  mk2(i,L) + dt * fmk2(i) + dt2 *
     X   (dfmk2_dmk2(i) * fmk2(i) + dfmk2_dv(i) * fv(i))
        hk2(i,L) =  hk2(i,L) + dt * fhk2(i) + dt2 *
     X   (dfhk2_dhk2(i) * fhk2(i) + dfhk2_dv(i) * fv(i))
        mkm(i,L) =  mkm(i,L) + dt * fmkm(i) + dt2 *
     X   (dfmkm_dmkm(i) * fmkm(i) + dfmkm_dv(i) * fv(i))
        mkc(i,L) =  mkc(i,L) + dt * fmkc(i) + dt2 *
     X   (dfmkc_dmkc(i) * fmkc(i) + dfmkc_dv(i) * fv(i))
        mkahp(i,L) = mkahp(i,L) + dt * fmkahp(i) + dt2 *
     X (dfmkahp_dmkahp(i)*fmkahp(i) + dfmkahp_dchi(i)*fchi(i))
        mcat(i,L) =  mcat(i,L) + dt * fmcat(i) + dt2 *
     X   (dfmcat_dmcat(i) * fmcat(i) + dfmcat_dv(i) * fv(i))
        hcat(i,L) =  hcat(i,L) + dt * fhcat(i) + dt2 *
     X   (dfhcat_dhcat(i) * fhcat(i) + dfhcat_dv(i) * fv(i))
        mcaL(i,L) =  mcaL(i,L) + dt * fmcaL(i) + dt2 *
     X   (dfmcaL_dmcaL(i) * fmcaL(i) + dfmcaL_dv(i) * fv(i))
        mar(i,L) =   mar(i,L) + dt * fmar(i) + dt2 *
     X   (dfmar_dmar(i) * fmar(i) + dfmar_dv(i) * fv(i))
c            endif
         end do

! Add membrane currents into membcurr for appropriate compartments
          do i = 1, 9
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do
          do i = 14, 21
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do
          do i = 26, 33
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do
          do i = 39, 68
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do

            end do
c Finish loop L = 1 to numcell

         field_sup = 0.d0
         field_deep = 0.d0

         do i = 1, 12
        field_sup = field_sup + membcurr(i) / dabs(100.d0 - depth(i))
        field_deep = field_deep + membcurr(i) / dabs(500.d0 - depth(i))
         end do

2001          CONTINUE

6000    END



C  SETS UP TABLES FOR RATE FUNCTIONS
       SUBROUTINE SCORT_SETUP_L2pyr   
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)
      INTEGER I,J,K
      real*8 minf, hinf, taum, tauh, V, Z, shift_hnaf,
     X  shift_mkdr,
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
C FOR VOLTAGE, RANGE IS -120 TO +40 MV (absol.), 0.25 MV RESOLUTION


       DO 1, I = 0, 640
          V = dble(I)
          V = (V / 4.d0) - 120.d0

c gNa
           minf = 1.d0/(1.d0 + dexp((-V-38.d0)/10.d0))
           if (v.le.-30.d0) then
            taum = .025d0 + .14d0*dexp((v+30.d0)/10.d0)
           else
            taum = .02d0 + .145d0*dexp((-v-30.d0)/10.d0)
           endif
c from principal c. data, Martina & Jonas 1997, tau x 0.5
c Note that minf about the same for interneuron & princ. cell.
           alpham_naf(i) = minf / taum
           betam_naf(i) = 1.d0/taum - alpham_naf(i)

            shift_hnaf =  0.d0
        hinf = 1.d0/(1.d0 +
     x     dexp((v + shift_hnaf + 62.9d0)/10.7d0))
        tauh = 0.15d0 + 1.15d0/(1.d0+dexp((v+37.d0)/15.d0))
c from princ. cell data, Martina & Jonas 1997, tau x 0.5
            alphah_naf(i) = hinf / tauh
            betah_naf(i) = 1.d0/tauh - alphah_naf(i)

          shift_mkdr = 0.d0
c delayed rectifier, non-inactivating
       minf = 1.d0/(1.d0+dexp((-v-shift_mkdr-29.5d0)/10.0d0))
            if (v.le.-10.d0) then
             taum = .25d0 + 4.35d0*dexp((v+10.d0)/10.d0)
            else
             taum = .25d0 + 4.35d0*dexp((-v-10.d0)/10.d0)
            endif
              alpham_kdr(i) = minf / taum
              betam_kdr(i) = 1.d0 /taum - alpham_kdr(i)
c from Martina, Schultz et al., 1998. See espec. Table 1.

c A current: Huguenard & McCormick 1992, J Neurophysiol (TCR)
            minf = 1.d0/(1.d0 + dexp((-v-60.d0)/8.5d0))
            hinf = 1.d0/(1.d0 + dexp((v+78.d0)/6.d0))
        taum = .185d0 + .5d0/(dexp((v+35.8d0)/19.7d0) +
     x                            dexp((-v-79.7d0)/12.7d0))
        if (v.le.-63.d0) then
         tauh = .5d0/(dexp((v+46.d0)/5.d0) +
     x                  dexp((-v-238.d0)/37.5d0))
        else
         tauh = 9.5d0
        endif
           alpham_ka(i) = minf/taum
           betam_ka(i) = 1.d0 / taum - alpham_ka(i)
           alphah_ka(i) = hinf / tauh
           betah_ka(i) = 1.d0 / tauh - alphah_ka(i)

c h-current (anomalous rectifier), Huguenard & McCormick, 1992
           minf = 1.d0/(1.d0 + dexp((v+75.d0)/5.5d0))
           taum = 1.d0/(dexp(-14.6d0 -0.086d0*v) +
     x                   dexp(-1.87 + 0.07d0*v))
           alpham_ar(i) = minf / taum
           betam_ar(i) = 1.d0 / taum - alpham_ar(i)

c K2 K-current, McCormick & Huguenard
             minf = 1.d0/(1.d0 + dexp((-v-10.d0)/17.d0))
             hinf = 1.d0/(1.d0 + dexp((v+58.d0)/10.6d0))
            taum = 4.95d0 + 0.5d0/(dexp((v-81.d0)/25.6d0) +
     x                  dexp((-v-132.d0)/18.d0))
            tauh = 60.d0 + 0.5d0/(dexp((v-1.33d0)/200.d0) +
     x                  dexp((-v-130.d0)/7.1d0))
             alpham_k2(i) = minf / taum
             betam_k2(i) = 1.d0/taum - alpham_k2(i)
             alphah_k2(i) = hinf / tauh
             betah_k2(i) = 1.d0 / tauh - alphah_k2(i)

c voltage part of C-current, using 1994 kinetics, shift 60 mV
              if (v.le.-10.d0) then
       alpham_kc(i) = (2.d0/37.95d0)*dexp((v+50.d0)/11.d0 -
     x                                     (v+53.5)/27.d0)
       betam_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)-alpham_kc(i)
               else
       alpham_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)
       betam_kc(i) = 0.d0
               endif

c high-threshold gCa, from 1994, with 60 mV shift & no inactivn.
            alpham_cal(i) = 1.6d0/(1.d0+dexp(-.072d0*(v-5.d0)))
            betam_cal(i) = 0.1d0 * ((v+8.9d0)/5.d0) /
     x          (dexp((v+8.9d0)/5.d0) - 1.d0)

c M-current, from plast.f, with 60 mV shift
        alpham_km(i) = .02d0/(1.d0+dexp((-v-20.d0)/5.d0))
        betam_km(i) = .01d0 * dexp((-v-43.d0)/18.d0)

c T-current, from Destexhe, Neubig et al., 1998
         minf = 1.d0/(1.d0 + dexp((-v-56.d0)/6.2d0))
         hinf = 1.d0/(1.d0 + dexp((v+80.d0)/4.d0))
         taum = 0.204d0 + .333d0/(dexp((v+15.8d0)/18.2d0) +
     x                  dexp((-v-131.d0)/16.7d0))
          if (v.le.-81.d0) then
         tauh = 0.333 * dexp((v+466.d0)/66.6d0)
          else
         tauh = 9.32d0 + 0.333d0*dexp((-v-21.d0)/10.5d0)
          endif
              alpham_cat(i) = minf / taum
              betam_cat(i) = 1.d0/taum - alpham_cat(i)
              alphah_cat(i) = hinf / tauh
              betah_cat(i) = 1.d0 / tauh - alphah_cat(i)

1        CONTINUE

         do  i = 0, 639

      dalpham_naf(i) = (alpham_naf(i+1)-alpham_naf(i))/.25d0
      dbetam_naf(i) = (betam_naf(i+1)-betam_naf(i))/.25d0
      dalphah_naf(i) = (alphah_naf(i+1)-alphah_naf(i))/.25d0
      dbetah_naf(i) = (betah_naf(i+1)-betah_naf(i))/.25d0
      dalpham_kdr(i) = (alpham_kdr(i+1)-alpham_kdr(i))/.25d0
      dbetam_kdr(i) = (betam_kdr(i+1)-betam_kdr(i))/.25d0
      dalpham_ka(i) = (alpham_ka(i+1)-alpham_ka(i))/.25d0
      dbetam_ka(i) = (betam_ka(i+1)-betam_ka(i))/.25d0
      dalphah_ka(i) = (alphah_ka(i+1)-alphah_ka(i))/.25d0
      dbetah_ka(i) = (betah_ka(i+1)-betah_ka(i))/.25d0
      dalpham_k2(i) = (alpham_k2(i+1)-alpham_k2(i))/.25d0
      dbetam_k2(i) = (betam_k2(i+1)-betam_k2(i))/.25d0
      dalphah_k2(i) = (alphah_k2(i+1)-alphah_k2(i))/.25d0
      dbetah_k2(i) = (betah_k2(i+1)-betah_k2(i))/.25d0
      dalpham_km(i) = (alpham_km(i+1)-alpham_km(i))/.25d0
      dbetam_km(i) = (betam_km(i+1)-betam_km(i))/.25d0
      dalpham_kc(i) = (alpham_kc(i+1)-alpham_kc(i))/.25d0
      dbetam_kc(i) = (betam_kc(i+1)-betam_kc(i))/.25d0
      dalpham_cat(i) = (alpham_cat(i+1)-alpham_cat(i))/.25d0
      dbetam_cat(i) = (betam_cat(i+1)-betam_cat(i))/.25d0
      dalphah_cat(i) = (alphah_cat(i+1)-alphah_cat(i))/.25d0
      dbetah_cat(i) = (betah_cat(i+1)-betah_cat(i))/.25d0
      dalpham_caL(i) = (alpham_cal(i+1)-alpham_cal(i))/.25d0
      dbetam_caL(i) = (betam_cal(i+1)-betam_cal(i))/.25d0
      dalpham_ar(i) = (alpham_ar(i+1)-alpham_ar(i))/.25d0
      dbetam_ar(i) = (betam_ar(i+1)-betam_ar(i))/.25d0
       end do
2      CONTINUE

         do i = 640, 640
      dalpham_naf(i) =  dalpham_naf(i-1)
      dbetam_naf(i) =  dbetam_naf(i-1)
      dalphah_naf(i) = dalphah_naf(i-1)
      dbetah_naf(i) = dbetah_naf(i-1)
      dalpham_kdr(i) =  dalpham_kdr(i-1)
      dbetam_kdr(i) =  dbetam_kdr(i-1)
      dalpham_ka(i) =  dalpham_ka(i-1)
      dbetam_ka(i) =  dbetam_ka(i-1)
      dalphah_ka(i) =  dalphah_ka(i-1)
      dbetah_ka(i) =  dbetah_ka(i-1)
      dalpham_k2(i) =  dalpham_k2(i-1)
      dbetam_k2(i) =  dbetam_k2(i-1)
      dalphah_k2(i) =  dalphah_k2(i-1)
      dbetah_k2(i) =  dbetah_k2(i-1)
      dalpham_km(i) =  dalpham_km(i-1)
      dbetam_km(i) =  dbetam_km(i-1)
      dalpham_kc(i) =  dalpham_kc(i-1)
      dbetam_kc(i) =  dbetam_kc(i-1)
      dalpham_cat(i) =  dalpham_cat(i-1)
      dbetam_cat(i) =  dbetam_cat(i-1)
      dalphah_cat(i) =  dalphah_cat(i-1)
      dbetah_cat(i) =  dbetah_cat(i-1)
      dalpham_caL(i) =  dalpham_caL(i-1)
      dbetam_caL(i) =  dbetam_caL(i-1)
      dalpham_ar(i) =  dalpham_ar(i-1)
      dbetam_ar(i) =  dbetam_ar(i-1)
       end do   

4000   END

        SUBROUTINE SCORTMAJ_L2pyr   
C BRANCHED ACTIVE DENDRITES
     X             (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM,depth,level)
c Conductances: leak gL, coupling g, delayed rectifier gKDR, A gKA,
c C gKC, AHP gKAHP, K2 gK2, M gKM, low thresh Ca gCAT, high thresh
c gCAL, fast Na gNAF, persistent Na gNAP, h or anom. rectif. gAR.
c Note VAR = equil. potential for anomalous rectifier.
c Soma = comp. 1; 10 dendrites each with 13 compartments, 6-comp. axon
c Drop "glc"-like terms, just using "gl"-like
c CAFOR corresponds to "phi" in Traub et al., 1994
c Consistent set of units: nF, mV, ms, nA, microS

       INTEGER, PARAMETER:: numcomp = 74
! numcomp here must be compatible with numcomp_suppyrRS in calling prog.
        REAL*8 C(numcomp),GL(numcomp), GAM(0:numcomp, 0:numcomp)
        REAL*8 GNAF(numcomp),GCAT(numcomp), GKAHP(numcomp)
        REAL*8 GKDR(numcomp),GKA(numcomp),GKC(numcomp)
        REAL*8 GK2(numcomp),GNAP(numcomp),GAR(numcomp)
        REAL*8 GKM(numcomp), gcal(numcomp), CDENS
        REAL*8 JACOB(numcomp,numcomp),RI_SD,RI_AXON,RM_SD,RM_AXON
        INTEGER LEVEL(numcomp)
        REAL*8 GNAF_DENS(0:12), GCAT_DENS(0:12), GKDR_DENS(0:12)
        REAL*8 GKA_DENS(0:12), GKC_DENS(0:12), GKAHP_DENS(0:12)
        REAL*8 GCAL_DENS(0:12), GK2_DENS(0:12), GKM_DENS(0:12)
        REAL*8 GNAP_DENS(0:12), GAR_DENS(0:12)
        REAL*8 RES, RINPUT, Z, ELEN(numcomp)
        REAL*8 RSOMA, PI, BETCHI(numcomp), CAFOR(numcomp)
        REAL*8 RAD(numcomp), LEN(numcomp), GAM1, GAM2
        REAL*8 RIN, D(numcomp), AREA(numcomp), RI
        INTEGER NEIGH(numcomp,10), NNUM(numcomp), i, j, k, it
C FOR ESTABLISHING TOPOLOGY OF COMPARTMENTS
        real*8 depth(12) ! depth in microns of levels 1-12, assuming soma 
! at depth 500 microns 

        depth(1) = 500.d0
        depth(2) = 550.d0
        depth(3) = 600.d0
        depth(4) = 650.d0
        depth(5) = 450.d0
        depth(6) = 400.d0
        depth(7) = 350.d0
        depth(8) = 300.d0
        depth(9) = 250.d0
        depth(10) = 200.d0
        depth(11) = 100.d0
        depth(12) =  50.d0

        RI_SD = 250.d0
        RM_SD = 50000.d0
        RI_AXON = 100.d0
        RM_AXON = 1000.d0
        CDENS = 0.9d0

        PI = 3.14159d0

       do i = 0, 12
        gnaf_dens(i) = 10.d0
       end do
c       gnaf_dens(0) = 400.d0
!       gnaf_dens(0) = 120.d0
        gnaf_dens(0) = 200.d0
        gnaf_dens(1) = 120.d0
        gnaf_dens(2) =  75.d0
        gnaf_dens(5) = 100.d0
        gnaf_dens(6) =  75.d0

       do i = 0, 12
        gkdr_dens(i) = 0.d0
       end do
c       gkdr_dens(0) = 400.d0
c       gkdr_dens(0) = 100.d0
c       gkdr_dens(0) = 170.d0
        gkdr_dens(0) = 250.d0
c       gkdr_dens(1) = 100.d0
        gkdr_dens(1) = 150.d0
        gkdr_dens(2) =  75.d0
        gkdr_dens(5) = 100.d0
        gkdr_dens(6) =  75.d0

        gnap_dens(0) = 0.d0
        do i = 1, 12
c         gnap_dens(i) = 0.0040d0 * gnaf_dens(i)
          gnap_dens(i) = 0.0010d0 * gnaf_dens(i)
c         gnap_dens(i) = 0.002d0 * gnaf_dens(i)
c         gnap_dens(i) = 0.0030d0 * gnaf_dens(i)
        end do

        gcat_dens(0) = 0.d0
        do i = 1, 12
c         gcat_dens(i) = 0.5d0
          gcat_dens(i) = 0.1d0
        end do

        gcaL_dens(0) = 0.d0
        do i = 1, 6
          gcaL_dens(i) = 0.5d0
c         gcaL_dens(i) = 0.060d0
        end do
        do i = 7, 12
          gcaL_dens(i) = 0.5d0
c         gcaL_dens(i) = 0.060d0
        end do

       do i = 0, 12
        gka_dens(i) = 2.d0
       end do
        gka_dens(0) =100.d0 ! NOTE
        gka_dens(1) = 30.d0
        gka_dens(5) = 30.d0

      do i = 0, 12
c        gkc_dens(i)  = 12.00d0
         gkc_dens(i)  =  0.00d0
c        gkc_dens(i)  =  2.00d0
c        gkc_dens(i)  =  7.00d0
      end do
         gkc_dens(0) =  0.00d0
c        gkc_dens(1) = 7.5d0
c        gkc_dens(1) = 12.d0
         gkc_dens(1) = 15.d0
c        gkc_dens(2) = 7.5d0
         gkc_dens(2) = 10.d0
         gkc_dens(5) = 7.5d0
         gkc_dens(6) = 7.5d0

c       gkm_dens(0) = 2.d0 ! 9 Nov. 2005, see scort-pan.f of today
        gkm_dens(0) = 8.d0 ! 9 Nov. 2005, see scort-pan.f of today
! Above suppresses doublets, but still allows FRB with appropriate
! gNaP, gKC, and rel_axonshift (e.g. 6 mV)
        do i = 1, 12
         gkm_dens(i) = 2.5d0 * 1.50d0
        end do

        do i = 0, 12
c       gk2_dens(i) = 1.d0
        gk2_dens(i) = 0.1d0
        end do
        gk2_dens(0) = 0.d0

        gkahp_dens(0) = 0.d0
        do i = 1, 12
c        gkahp_dens(i) = 0.200d0
         gkahp_dens(i) = 0.100d0
c        gkahp_dens(i) = 0.050d0
        end do

        gar_dens(0) = 0.d0
        do i = 1, 12
         gar_dens(i) = 0.25d0
        end do

c       WRITE   (6,9988)
9988    FORMAT(2X,'I',4X,'NADENS',' CADENS(T)',' KDRDEN',' KAHPDE',
     X     ' KCDENS',' KADENS')
        DO 9989, I = 0, 12
c         WRITE (6,9990) I, gnaf_dens(i), gcat_dens(i), gkdr_dens(i),
c    X  gkahp_dens(i), gkc_dens(i), gka_dens(i)
9990    FORMAT(2X,I2,2X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2)
9989    CONTINUE


        level(1) = 1
        do i = 2, 13
         level(i) = 2
        end do
        do i = 14, 25
           level(i) = 3
        end do
        do i = 26, 37
           level(i) = 4
        end do
        level(38) = 5
        level(39) = 6
        level(40) = 7
        level(41) = 8
        level(42) = 8
        level(43) = 9
        level(44) = 9
        do i = 45, 52
           level(i) = 10
        end do
        do i = 53, 60
           level(i) = 11
        end do
        do i = 61, 68
           level(i) = 12
        end do

        do i =  69, 74
         level(i) = 0
        end do

c connectivity of axon
        nnum( 69) = 2
        nnum( 70) = 3
        nnum( 71) = 3
        nnum( 73) = 3
        nnum( 72) = 1
        nnum( 74) = 1
         neigh(69,1) =  1
         neigh(69,2) = 70
         neigh(70,1) = 69
         neigh(70,2) = 71
         neigh(70,3) = 73
         neigh(71,1) = 70
         neigh(71,2) = 72
         neigh(71,3) = 73
         neigh(73,1) = 70
         neigh(73,2) = 71
         neigh(73,3) = 74
         neigh(72,1) = 71
         neigh(74,1) = 73

c connectivity of SD part
          nnum(1) = 10
          neigh(1,1) = 69
          neigh(1,2) =  2
          neigh(1,3) =  3
          neigh(1,4) =  4
          neigh(1,5) =  5
          neigh(1,6) =  6
          neigh(1,7) =  7
          neigh(1,8) =  8
          neigh(1,9) =  9
          neigh(1,10) = 38

          do i = 2, 9
           nnum(i) = 2
           neigh(i,1) = 1
           neigh(i,2) = i + 12
          end do

          do i = 14, 21
            nnum(i) = 2
            neigh(i,1) = i - 12
            neigh(i,2) = i + 12
          end do

          do i = 26, 33
            nnum(i) = 1
            neigh(i,1) = i - 12
          end do

          do i = 10, 13
            nnum(i) = 2
            neigh(i,1) = 38
            neigh(i,2) = i + 12
          end do

          do i = 22, 25
            nnum(i) = 2
            neigh(i,1) = i - 12
            neigh(i,2) = i + 12
          end do

          do i = 34, 37
            nnum(i) = 1
            neigh(i,1) = i - 12
          end do

          nnum(38) = 6
          neigh(38,1) = 1
          neigh(38,2) = 39
          neigh(38,3) = 10
          neigh(38,4) = 11
          neigh(38,5) = 12
          neigh(38,6) = 13

          nnum(39) = 2
          neigh(39,1) = 38
          neigh(39,2) = 40

          nnum(40) = 3
          neigh(40,1) = 39
          neigh(40,2) = 41
          neigh(40,3) = 42

          nnum(41) = 3
          neigh(41,1) = 40
          neigh(41,2) = 42
          neigh(41,3) = 43

          nnum(42) = 3
          neigh(42,1) = 40
          neigh(42,2) = 41
          neigh(42,3) = 44

           nnum(43) = 5
           neigh(43,1) = 41
           neigh(43,2) = 45
           neigh(43,3) = 46
           neigh(43,4) = 47
           neigh(43,5) = 48

           nnum(44) = 5
           neigh(44,1) = 42
           neigh(44,2) = 49
           neigh(44,3) = 50
           neigh(44,4) = 51
           neigh(44,5) = 52

           nnum(45) = 5
           neigh(45,1) = 43
           neigh(45,2) = 53
           neigh(45,3) = 46
           neigh(45,4) = 47
           neigh(45,5) = 48

           nnum(46) = 5
           neigh(46,1) = 43
           neigh(46,2) = 54
           neigh(46,3) = 45
           neigh(46,4) = 47
           neigh(46,5) = 48

           nnum(47) = 5
           neigh(47,1) = 43
           neigh(47,2) = 55
           neigh(47,3) = 45
           neigh(47,4) = 46
           neigh(47,5) = 48

           nnum(48) = 5
           neigh(48,1) = 43
           neigh(48,2) = 56
           neigh(48,3) = 45
           neigh(48,4) = 46
           neigh(48,5) = 47

           nnum(49) = 5
           neigh(49,1) = 44
           neigh(49,2) = 57
           neigh(49,3) = 50
           neigh(49,4) = 51
           neigh(49,5) = 52

           nnum(50) = 5
           neigh(50,1) = 44
           neigh(50,2) = 58
           neigh(50,3) = 49
           neigh(50,4) = 51
           neigh(50,5) = 52

           nnum(51) = 5
           neigh(51,1) = 44
           neigh(51,2) = 59
           neigh(51,3) = 49
           neigh(51,4) = 50
           neigh(51,5) = 52

           nnum(52) = 5
           neigh(52,1) = 44
           neigh(52,2) = 60
           neigh(52,3) = 49
           neigh(52,4) = 51
           neigh(52,5) = 50

          do i = 53, 60
           nnum(i) = 2
           neigh(i,1) = i - 8
           neigh(i,2) = i + 8
          end do

          do i = 61, 68
           nnum(i) = 1
           neigh(i,1) = i - 8
          end do

c        DO 332, I = 1, 74
         DO I = 1, 74
c          WRITE(6,3330) I, NEIGH(I,1),NEIGH(I,2),NEIGH(I,3),NEIGH(I,4),
c    X NEIGH(I,5),NEIGH(I,6),NEIGH(I,7),NEIGH(I,8),NEIGH(I,9),
c    X NEIGH(I,10)
3330     FORMAT(2X,11I5)
         END DO
332      CONTINUE
c         DO 858, I = 1, 74
          DO I = 1, 74
c          DO 858, J = 1, NNUM(I)
           DO J = 1, NNUM(I)
            K = NEIGH(I,J)
            IT = 0
c           DO 859, L = 1, NNUM(K)
            DO  L = 1, NNUM(K)
             IF (NEIGH(K,L).EQ.I) IT = 1
            END DO
859         CONTINUE
             IF (IT.EQ.0) THEN
c             WRITE(6,8591) I, K
8591          FORMAT(' ASYMMETRY IN NEIGH MATRIX ',I4,I4)
              STOP
             ENDIF
          END DO
          END DO
858       CONTINUE

c length and radius of axonal compartments
c Note shortened "initial segment"
          len(69) = 25.d0
          do i = 70, 74
            len(i) = 50.d0
          end do
          rad( 69) = 0.90d0
c         rad( 69) = 0.80d0
          rad( 70) = 0.7d0
          do i = 71, 74
           rad(i) = 0.5d0
          end do

c  length and radius of SD compartments
          len(1) = 15.d0
          rad(1) =  8.d0

          do i = 2, 68
           len(i) = 50.d0
          end do

          do i = 2, 37
            rad(i) = 0.5d0
          end do

          z = 4.0d0
          rad(38) = z
          rad(39) = 0.9d0 * z
          rad(40) = 0.8d0 * z
          rad(41) = 0.5d0 * z
          rad(42) = 0.5d0 * z
          rad(43) = 0.5d0 * z
          rad(44) = 0.5d0 * z
          do i = 45, 68
           rad(i) = 0.2d0 * z
          end do


c       WRITE(6,919)
919     FORMAT('COMPART.',' LEVEL ',' RADIUS ',' LENGTH(MU)')
c       DO 920, I = 1, 74
c920      WRITE(6,921) I, LEVEL(I), RAD(I), LEN(I)
921     FORMAT(I3,5X,I2,3X,F6.2,1X,F6.1,2X,F4.3)

        DO 120, I = 1, 74
          AREA(I) = 2.d0 * PI * RAD(I) * LEN(I)
      if((i.gt.1).and.(i.le.68)) area(i) = 2.d0 * area(i)
C    CORRECTION FOR CONTRIBUTION OF SPINES TO AREA
          K = LEVEL(I)
          C(I) = CDENS * AREA(I) * (1.D-8)

           if (k.ge.1) then
          GL(I) = (1.D-2) * AREA(I) / RM_SD
           else
          GL(I) = (1.D-2) * AREA(I) / RM_AXON
           endif

          GNAF(I) = GNAF_DENS(K) * AREA(I) * (1.D-5)
          GNAP(I) = GNAP_DENS(K) * AREA(I) * (1.D-5)
          GCAT(I) = GCAT_DENS(K) * AREA(I) * (1.D-5)
          GKDR(I) = GKDR_DENS(K) * AREA(I) * (1.D-5)
          GKA(I) = GKA_DENS(K) * AREA(I) * (1.D-5)
          GKC(I) = GKC_DENS(K) * AREA(I) * (1.D-5)
          GKAHP(I) = GKAHP_DENS(K) * AREA(I) * (1.D-5)
          GKM(I) = GKM_DENS(K) * AREA(I) * (1.D-5)
          GCAL(I) = GCAL_DENS(K) * AREA(I) * (1.D-5)
          GK2(I) = GK2_DENS(K) * AREA(I) * (1.D-5)
          GAR(I) = GAR_DENS(K) * AREA(I) * (1.D-5)
c above conductances should be in microS
120           continue

         Z = 0.d0
c        DO 1019, I = 2, 68
         DO I = 2, 68
           Z = Z + AREA(I)
         END DO
1019     CONTINUE
c        WRITE(6,1020) Z
1020     FORMAT(2X,' TOTAL DENDRITIC AREA ',F7.0)

c       DO 140, I = 1, 74
        DO I = 1, 74
c       DO 140, K = 1, NNUM(I)
        DO K = 1, NNUM(I)
         J = NEIGH(I,K)
           if (level(i).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM1 =100.d0 * PI * RAD(I) * RAD(I) / ( RI * LEN(I) )

           if (level(j).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM2 =100.d0 * PI * RAD(J) * RAD(J) / ( RI * LEN(J) )
         GAM(I,J) = 2.d0/( (1.d0/GAM1) + (1.d0/GAM2) )
	 END DO
	 END DO

140     CONTINUE
c gam computed in microS

c       DO 299, I = 1, 74
        DO I = 1, 74
299       BETCHI(I) = .05d0
        END DO
        BETCHI( 1) =  .01d0

c       DO 300, I = 1, 74
        DO I = 1, 74
c300     D(I) = 2.D-4
300     D(I) = 5.D-4
        END DO
c       DO 301, I = 1, 74
        DO I = 1, 74
         IF (LEVEL(I).EQ.1) D(I) = 2.D-3
        END DO
301     CONTINUE
C  NOTE NOTE NOTE  (DIFFERENT FROM SWONG)


c      DO 160, I = 1, 74
       DO I = 1, 74
160     CAFOR(I) = 5200.d0 / (AREA(I) * D(I))
       END DO
C     NOTE CORRECTION

c       do 200, i = 1, 74
        do i = 1, numcomp
200     C(I) = 1000.d0 * C(I)
        end do
C     TO GO FROM MICROF TO NF.

c     DO 909, I = 1, 74
      DO I = 1, numcomp
       JACOB(I,I) = - GL(I)
c     DO 909, J = 1, NNUM(I)
      DO J = 1, NNUM(I)
         K = NEIGH(I,J)
         IF (I.EQ.K) THEN
c            WRITE(6,510) I
510          FORMAT(' UNEXPECTED SYMMETRY IN NEIGH ',I4)
         ENDIF
         JACOB(I,K) = GAM(I,K)
         JACOB(I,I) = JACOB(I,I) - GAM(I,K)
       END DO
       END DO
909   CONTINUE

c 15 Jan. 2001: make correction for c(i)
          do i = 1, numcomp
          do j = 1, numcomp
             jacob(i,j) = jacob(i,j) / c(i)
          end do
          end do

c      DO 500, I = 1, 74
       DO I = 1, 74
c       WRITE (6,501) I,C(I)
501     FORMAT(1X,I3,' C(I) = ',F7.4)
       END DO
500     CONTINUE
        END

c 22 Aug 2019, start with suppyrRS integration subroutine from
c son_of_groucho, and use for L3pyr in piriform simulations.
c Need to change field variables and depth definitions,
c  and perhaps alter compartment dimensions.
c 11 Sept 2006, start with /interact/integrate_suppyrRSXP.f & add GABA-B
! 7 Nov. 2005: modify integrate_suppyrRSX.f to allow for Colbert-Pan axon.
!29 July 2005: modify groucho/integrate_suppyrRS.f, for a separate
! call for initialization, and to integrate only selected cells.
! Integration routine for suppyrRS cells
! Routine adapted from scortn in supergj.f
c      SUBROUTINE INTEGRATE_suppyrRSXPB (O, time, numcell,     
       SUBROUTINE INTEGRATE_L3pyr       (O, time, numcell,     
     &    V, curr, initialize, firstcell, lastcell,
     & gAMPA, gNMDA, gGABA_A, gGABA_B,
     & Mg, 
     & gapcon  ,totaxgj   ,gjtable, dt,
     &  chi,mnaf,mnap,
     &  hnaf,mkdr,mka,
     &  hka,mk2,hk2,
     &  mkm,mkc,mkahp,
     &  mcat,hcat,mcal,
     &  mar,field_sup,field_deep,rel_axonshift)

       SAVE

       INTEGER, PARAMETER:: numcomp = 74
! numcomp here must be compatible with numcomp_suppyrRS in calling prog.
       INTEGER  numcell, num_other
       INTEGER initialize, firstcell, lastcell
       INTEGER J1, I, J, K, K1, K2, K3, L, L1, O
       REAL*8 c(numcomp), curr(numcomp,numcell)
       REAL*8  Z, Z1, Z2, Z3, Z4, DT, time
       integer totaxgj, gjtable(totaxgj,4)
       real*8 gapcon, gAMPA(numcomp,numcell),
     &        gNMDA(numcomp,numcell), gGABA_A(numcomp,numcell),
     &        gGABA_B(numcomp,numcell)
       real*8 Mg, V(numcomp,numcell), rel_axonshift

c CINV is 1/C, i.e. inverse capacitance
       real*8 chi(numcomp,numcell),
     & mnaf(numcomp,numcell),mnap(numcomp,numcell),
     x hnaf(numcomp,numcell), mkdr(numcomp,numcell),
     x mka(numcomp,numcell),hka(numcomp,numcell),
     x mk2(numcomp,numcell), cinv(numcomp),
     x hk2(numcomp,numcell),mkm(numcomp,numcell),
     x mkc(numcomp,numcell),mkahp(numcomp,numcell),
     x mcat(numcomp,numcell),hcat(numcomp,numcell),
     x mcal(numcomp,numcell), betchi(numcomp),
     x mar(numcomp,numcell),jacob(numcomp,numcomp),
     x gam(0: numcomp,0: numcomp),gL(numcomp),gnaf(numcomp),
     x gnap(numcomp),gkdr(numcomp),gka(numcomp),
     x gk2(numcomp),gkm(numcomp),
     x gkc(numcomp),gkahp(numcomp),
     x gcat(numcomp),gcaL(numcomp),gar(numcomp),
     x cafor(numcomp)
       real*8
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
       real*8 vL(numcomp),vk(numcomp),vna,var,vca,vgaba_a
       real*8 depth(12), membcurr(12), field_sup, field_deep
       integer level(numcomp)

        INTEGER NEIGH(numcomp,10), NNUM(numcomp)
        INTEGER igap1, igap2
c the f's are the functions giving 1st derivatives for evolution of
c the differential equations for the voltages (v), calcium (chi), and
c other state variables.
       real*8 fv(numcomp), fchi(numcomp),
     x fmnaf(numcomp),fhnaf(numcomp),fmkdr(numcomp),
     x fmka(numcomp),fhka(numcomp),fmk2(numcomp),
     x fhk2(numcomp),fmnap(numcomp),
     x fmkm(numcomp),fmkc(numcomp),fmkahp(numcomp),
     x fmcat(numcomp),fhcat(numcomp),fmcal(numcomp),
     x fmar(numcomp)

c below are for calculating the partial derivatives
       real*8 dfv_dv(numcomp,numcomp), dfv_dchi(numcomp),
     x  dfv_dmnaf(numcomp),  dfv_dmnap(numcomp),
     x  dfv_dhnaf(numcomp),dfv_dmkdr(numcomp),
     x  dfv_dmka(numcomp),dfv_dhka(numcomp),
     x  dfv_dmk2(numcomp),dfv_dhk2(numcomp),
     x  dfv_dmkm(numcomp),dfv_dmkc(numcomp),
     x  dfv_dmkahp(numcomp),dfv_dmcat(numcomp),
     x  dfv_dhcat(numcomp),dfv_dmcal(numcomp),
     x  dfv_dmar(numcomp)

        real*8 dfchi_dv(numcomp), dfchi_dchi(numcomp),
     x dfmnaf_dmnaf(numcomp), dfmnaf_dv(numcomp),
     x dfhnaf_dhnaf(numcomp),
     x dfmnap_dmnap(numcomp), dfmnap_dv(numcomp),
     x dfhnaf_dv(numcomp),dfmkdr_dmkdr(numcomp),
     x dfmkdr_dv(numcomp),
     x dfmka_dmka(numcomp),dfmka_dv(numcomp),
     x dfhka_dhka(numcomp),dfhka_dv(numcomp),
     x dfmk2_dmk2(numcomp),dfmk2_dv(numcomp),
     x dfhk2_dhk2(numcomp),dfhk2_dv(numcomp),
     x dfmkm_dmkm(numcomp),dfmkm_dv(numcomp),
     x dfmkc_dmkc(numcomp),dfmkc_dv(numcomp),
     x dfmcat_dmcat(numcomp),dfmcat_dv(numcomp),dfhcat_dhcat(numcomp),
     x dfhcat_dv(numcomp),dfmcal_dmcal(numcomp),dfmcal_dv(numcomp),
     x dfmar_dmar(numcomp),dfmar_dv(numcomp),dfmkahp_dchi(numcomp),
     x dfmkahp_dmkahp(numcomp), dt2

       REAL*8 OPEN(numcomp),gamma(numcomp),gamma_prime(numcomp)
c gamma is function of chi used in calculating KC conductance
       REAL*8 alpham_ahp(numcomp), alpham_ahp_prime(numcomp)
       REAL*8 gna_tot(numcomp),gk_tot(numcomp),gca_tot(numcomp)
       REAL*8 gca_high(numcomp), gar_tot(numcomp)
c this will be gCa conductance corresponding to high-thresh channels

       real*8 persistentNa_shift, fastNa_shift_SD,
     x   fastNa_shift_axon

       REAL*8 A, BB1, BB2  ! params. for FNMDA.f


c          if (O.eq.1) then
           if (initialize.eq.0) then
c do initialization

c Program fnmda assumes A, BB1, BB2 defined in calling program
c as follows:
         A = DEXP(-2.847d0)
         BB1 = DEXP(-.693d0)
         BB2 = DEXP(-3.101d0)

c       goto 4000
       CALL   SCORT_SETUP_L3pyr   
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)

        CALL SCORTMAJ_L3pyr   
     X             (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM,depth,level)

          do i = 1, numcomp
             cinv(i) = 1.d0 / c(i)
          end do
4000      CONTINUE

           do i = 1, numcomp
          vL(i) = -70.d0
          vK(i) = -95.d0
           end do

        VNA = 50.d0
        VCA = 125.d0
        VAR = -43.d0
        VAR = -35.d0
c -43 mV from Huguenard & McCormick
        VGABA_A = -81.d0
c       write(6,901) VNa, VCa, VK(1), O
901     format('VNa =',f6.2,' VCa =',f6.2,' VK =',f6.2,
     &   ' O = ',i3)

c ? initialize membrane state variables?
         do L = 1, numcell  
         do i = 1, numcomp
        v(i,L) = VL(i)
	chi(i,L) = 0.d0
	mnaf(i,L) = 0.d0
	mkdr(i,L) = 0.d0
	mk2(i,L) = 0.d0
	mkm(i,L) = 0.d0
	mkc(i,L) = 0.d0
	mkahp(i,L) = 0.d0
	mcat(i,L) = 0.d0
	mcal(i,L) = 0.d0
         end do
         end do

          do L = 1, numcell
        k1 = idnint (4.d0 * (v(1,L) + 120.d0))

            do i = 1, numcomp
      hnaf(i,L) = alphah_naf(k1)/(alphah_naf(k1)
     &       +betah_naf(k1))
      hka(i,L) = alphah_ka(k1)/(alphah_ka(k1)
     &                               +betah_ka(k1))
      hk2(i,L) = alphah_k2(k1)/(alphah_k2(k1)
     &                                +betah_k2(k1))
      hcat(i,L)=alphah_cat(k1)/(alphah_cat(k1)
     &                                +betah_cat(k1))
c     mar=alpham_ar(k1)/(alpham_ar(k1)+betam_ar(k1))
      mar(i,L) = .25d0
             end do
           end do


             do i = 1, numcomp
	    open(i) = 0.d0
            gkm(i) = 2.d0 * gkm(i)
             end do

         do i = 1, 68
c          gnaf(i) = 0.8d0 * 1.25d0 * gnaf(i) ! factor of 0.8 added 19 Nov. 2005
c          gnaf(i) = 0.9d0 * 1.25d0 * gnaf(i) ! Back to 0.9, 29 Nov. 2005
           gnaf(i) = 0.6d0 * 1.25d0 * gnaf(i) ! 
! NOTE THAT THERE IS QUESTION OF HOW TO COMPARE BEHAVIOR OF PYRAMID IN NETWORK WITH
! SIMULATIONS OF SINGLE CELL.  IN FORMER CASE, THERE IS LARGE AXONAL SHUNT THROUGH
! gj(s), NOT PRESENT IN SINGLE CELL MODEL.  THEREFORE, HIGHER AXONAL gNa MIGHT BE
! NECESSARY FOR SPIKE PROPAGATION.
c          gnaf(i) = 0.9d0 * 1.25d0 * gnaf(i) ! factor of 0.9 added 20 Nov. 2005
           gkdr(i) = 1.25d0 * gkdr(i)
         end do
 
c Perhaps reduce fast gNa on IS
          gnaf(69) = 1.00d0 * gnaf(69)
c         gnaf(69) = 0.25d0 * gnaf(69)
          gnaf(70) = 1.00d0 * gnaf(70)
c         gnaf(70) = 0.25d0 * gnaf(70)

c Perhaps reduce coupling between soma and IS
c         gam(1,69) = 0.15d0 * gam(1,69)
c         gam(69,1) = 0.15d0 * gam(69,1)

               z1 = 0.0d0
c              z2 = 1.2d0 ! value 1.2 tried Feb. 21, 2013
               z2 = 1.5d0 ! value 1.2 tried Feb. 21, 2013
               z3 = 1.0d0
c              z3 = 0.0d0 ! Note reduction from 0.4, to prevent
c slow hyperpolarization that seems to mess up gamma.
               z4 = 0.3d0
c RS cell
             do i = 1, numcomp
              gnap(i) = z1 * gnap(i)
              gkc (i) = z2 * gkc (i)
              gkahp(i) = z3 * gkahp(i)
              gkm (i) = z4 * gkm (i)
             end do

              goto 6000

          endif
c End initialization

          do i = 1, 12
           membcurr(i) = 0.d0
          end do

c                  goto 2001


c             do L = 1, numcell
              do L = firstcell, lastcell

	  do i = 1, numcomp
	  do j = 1, nnum(i)
	   if (neigh(i,j).gt.numcomp) then
          write(6,433) i, j, L
433       format(' ls ',3x,3i5)
           endif
	end do
	end do

       DO I = 1, numcomp
          FV(I) = -GL(I) * (V(I,L) - VL(i)) * cinv(i)
          DO J = 1, NNUM(I)
             K = NEIGH(I,J)
302     FV(I) = FV(I) + GAM(I,K) * (V(K,L) - V(I,L)) * cinv(i)
           END DO
       END DO
301    CONTINUE


       CALL FNMDA (V, OPEN, numcell, numcomp, MG, L,
     &                 A, BB1, BB2)

      DO I = 1, numcomp
       FV(I) = FV(I) + ( CURR(I,L)
     X   - (gampa(I,L) + open(i) * gnmda(I,L))*V(I,L)
     X   - ggaba_a(I,L)*(V(I,L)-Vgaba_a) 
     X   - ggaba_b(I,L)*(V(I,L)-VK(i)  ) ) * cinv(i)
c above assumes equil. potential for AMPA & NMDA = 0 mV
      END DO
421      continue

       do m = 1, totaxgj
        if (gjtable(m,1).eq.L) then
         L1 = gjtable(m,3)
         igap1 = gjtable(m,2)
         igap2 = gjtable(m,4)
 	fv(igap1) = fv(igap1) + gapcon *
     &   (v(igap2,L1) - v(igap1,L)) * cinv(igap1)
        else if (gjtable(m,3).eq.L) then
         L1 = gjtable(m,1)
         igap1 = gjtable(m,4)
         igap2 = gjtable(m,2)
 	fv(igap1) = fv(igap1) + gapcon *
     &   (v(igap2,L1) - v(igap1,L)) * cinv(igap1)
        endif
       end do ! do m


       do i = 1, numcomp
        gamma(i) = dmin1 (1.d0, .004d0 * chi(i,L))
        if (chi(i,L).le.250.d0) then
          gamma_prime(i) = .004d0
        else
          gamma_prime(i) = 0.d0
        endif
c         endif
       end do

      DO I = 1, numcomp
       gna_tot(i) = gnaf(i) * (mnaf(i,L)**3) * hnaf(i,L) +
     x     gnap(i) * mnap(i,L)
       gk_tot(i) = gkdr(i) * (mkdr(i,L)**4) +
     x             gka(i)  * (mka(i,L)**4) * hka(i,L) +
     x             gk2(i)  * mk2(i,L) * hk2(i,L) +
     x             gkm(i)  * mkm(i,L) +
     x             gkc(i)  * mkc(i,L) * gamma(i) +
     x             gkahp(i)* mkahp(i,L)
       gca_tot(i) = gcat(i) * (mcat(i,L)**2) * hcat(i,L) +
     x              gcaL(i) * (mcaL(i,L)**2)
       gca_high(i) =
     x              gcaL(i) * (mcaL(i,L)**2)
       gar_tot(i) = gar(i) * mar(i,L)


       FV(I) = FV(I) - ( gna_tot(i) * (v(i,L) - vna)
     X  + gk_tot(i) * (v(i,L) - vK(i))
     X  + gca_tot(i) * (v(i,L) - vCa)
     X  + gar_tot(i) * (v(i,L) - var) ) * cinv(i)
c        endif
       END DO
88           continue

         do i = 1, numcomp
         do j = 1, numcomp
          if (i.ne.j) then
            dfv_dv(i,j) = jacob(i,j)
          else
            dfv_dv(i,j) = jacob(i,i) - cinv(i) *
     X  (gna_tot(i) + gk_tot(i) + gca_tot(i) + gar_tot(i)
     X   + ggaba_a(i,L) + ggaba_b(i,L) + gampa(i,L)
     X   + open(i) * gnmda(I,L) )
          endif
         end do
         end do

           do i = 1, numcomp
        dfv_dchi(i)  = - cinv(i) * gkc(i) * mkc(i,L) *
     x                     gamma_prime(i) * (v(i,L)-vK(i))
        dfv_dmnaf(i) = -3.d0 * cinv(i) * (mnaf(i,L)**2) *
     X    (gnaf(i) * hnaf(i,L)          ) * (v(i,L) - vna)
        dfv_dmnap(i) = - cinv(i) *
     X    (               gnap(i)) * (v(i,L) - vna)
        dfv_dhnaf(i) = - cinv(i) * gnaf(i) * (mnaf(i,L)**3) *
     X                    (v(i,L) - vna)
        dfv_dmkdr(i) = -4.d0 * cinv(i) * gkdr(i) * (mkdr(i,L)**3)
     X                   * (v(i,L) - vK(i))
        dfv_dmka(i)  = -4.d0 * cinv(i) * gka(i) * (mka(i,L)**3) *
     X                   hka(i,L) * (v(i,L) - vK(i))
        dfv_dhka(i)  = - cinv(i) * gka(i) * (mka(i,L)**4) *
     X                    (v(i,L) - vK(i))
      dfv_dmk2(i) = - cinv(i) * gk2(i) * hk2(i,L) * (v(i,L)-vK(i))
      dfv_dhk2(i) = - cinv(i) * gk2(i) * mk2(i,L) * (v(i,L)-vK(i))
      dfv_dmkm(i) = - cinv(i) * gkm(i) * (v(i,L) - vK(i))
      dfv_dmkc(i) = - cinv(i)*gkc(i) * gamma(i) * (v(i,L)-vK(i))
        dfv_dmkahp(i)= - cinv(i) * gkahp(i) * (v(i,L) - vK(i))
        dfv_dmcat(i)  = -2.d0 * cinv(i) * gcat(i) * mcat(i,L) *
     X                    hcat(i,L) * (v(i,L) - vCa)
        dfv_dhcat(i) = - cinv(i) * gcat(i) * (mcat(i,L)**2) *
     X                  (v(i,L) - vCa)
        dfv_dmcal(i) = -2.d0 * cinv(i) * gcal(i) * mcal(i,L) *
     X                      (v(i,L) - vCa)
        dfv_dmar(i) = - cinv(i) * gar(i) * (v(i,L) - var)
            end do

         do i = 1, numcomp
          fchi(i) = - cafor(i) * gca_high(i) * (v(i,L) - vca)
     x       - betchi(i) * chi(i,L)
          dfchi_dv(i) = - cafor(i) * gca_high(i)
          dfchi_dchi(i) = - betchi(i)
         end do

       do i = 1, numcomp
c Note possible increase in rate at which AHP current develops
c       alpham_ahp(i) = dmin1(0.2d-4 * chi(i,L),0.01d0)
        alpham_ahp(i) = dmin1(1.0d-4 * chi(i,L),0.01d0)
        if (chi(i,L).le.500.d0) then
c         alpham_ahp_prime(i) = 0.2d-4
          alpham_ahp_prime(i) = 1.0d-4
        else
          alpham_ahp_prime(i) = 0.d0
        endif
       end do

       do i = 1, numcomp
        fmkahp(i) = alpham_ahp(i) * (1.d0 - mkahp(i,L))
c    x                  -.001d0 * mkahp(i,L)
     x                  -.010d0 * mkahp(i,L)
c       dfmkahp_dmkahp(i) = - alpham_ahp(i) - .001d0
        dfmkahp_dmkahp(i) = - alpham_ahp(i) - .010d0
        dfmkahp_dchi(i) = alpham_ahp_prime(i) *
     x                     (1.d0 - mkahp(i,L))
       end do

          do i = 1, numcomp

       K1 = IDNINT ( 4.d0 * (V(I,L) + 120.d0) )
       IF (K1.GT.640) K1 = 640
       IF (K1.LT.  0) K1 =   0

c      persistentNa_shift =  0.d0
c      persistentNa_shift =  8.d0
       persistentNa_shift = 10.d0
       K2 = IDNINT ( 4.d0 * (V(I,L)+persistentNa_shift+ 120.d0) )
       IF (K2.GT.640) K2 = 640
       IF (K2.LT.  0) K2 =   0

c            fastNa_shift = -2.0d0
c            fastNa_shift = -2.5d0
             fastNa_shift_SD = -3.5d0
             fastNa_shift_axon = fastNa_shift_SD + rel_axonshift 
       K0 = IDNINT ( 4.d0 * (V(I,L)+  fastNa_shift_SD+ 120.d0) )
       IF (K0.GT.640) K0 = 640
       IF (K0.LT.  0) K0 =   0
       K3 = IDNINT ( 4.d0 * (V(I,L)+  fastNa_shift_axon+ 120.d0) )
       IF (K3.GT.640) K3 = 640
       IF (K3.LT.  0) K3 =   0

         if (i.le.68) then   ! FOR SD
        fmnaf(i) = alpham_naf(k0) * (1.d0 - mnaf(i,L)) -
     X              betam_naf(k0) * mnaf(i,L)
        fhnaf(i) = alphah_naf(k0) * (1.d0 - hnaf(i,L)) -
     X              betah_naf(k0) * hnaf(i,L)
         else  ! for axon
        fmnaf(i) = alpham_naf(k3) * (1.d0 - mnaf(i,L)) -
     X              betam_naf(k3) * mnaf(i,L)
        fhnaf(i) = alphah_naf(k3) * (1.d0 - hnaf(i,L)) -
     X              betah_naf(k3) * hnaf(i,L)
         endif
        fmnap(i) = alpham_naf(k2) * (1.d0 - mnap(i,L)) -
     X              betam_naf(k2) * mnap(i,L)
        fmkdr(i) = alpham_kdr(k1) * (1.d0 - mkdr(i,L)) -
     X              betam_kdr(k1) * mkdr(i,L)
        fmka(i)  = alpham_ka (k1) * (1.d0 - mka(i,L)) -
     X              betam_ka (k1) * mka(i,L)
        fhka(i)  = alphah_ka (k1) * (1.d0 - hka(i,L)) -
     X              betah_ka (k1) * hka(i,L)
        fmk2(i)  = alpham_k2 (k1) * (1.d0 - mk2(i,L)) -
     X              betam_k2 (k1) * mk2(i,L)
        fhk2(i)  = alphah_k2 (k1) * (1.d0 - hk2(i,L)) -
     X              betah_k2 (k1) * hk2(i,L)
        fmkm(i)  = alpham_km (k1) * (1.d0 - mkm(i,L)) -
     X              betam_km (k1) * mkm(i,L)
        fmkc(i)  = alpham_kc (k1) * (1.d0 - mkc(i,L)) -
     X              betam_kc (k1) * mkc(i,L)
        fmcat(i) = alpham_cat(k1) * (1.d0 - mcat(i,L)) -
     X              betam_cat(k1) * mcat(i,L)
        fhcat(i) = alphah_cat(k1) * (1.d0 - hcat(i,L)) -
     X              betah_cat(k1) * hcat(i,L)
        fmcaL(i) = alpham_caL(k1) * (1.d0 - mcaL(i,L)) -
     X              betam_caL(k1) * mcaL(i,L)
        fmar(i)  = alpham_ar (k1) * (1.d0 - mar(i,L)) -
     X              betam_ar (k1) * mar(i,L)

       dfmnaf_dv(i) = dalpham_naf(k0) * (1.d0 - mnaf(i,L)) -
     X                  dbetam_naf(k0) * mnaf(i,L)
       dfmnap_dv(i) = dalpham_naf(k2) * (1.d0 - mnap(i,L)) -
     X                  dbetam_naf(k2) * mnap(i,L)
       dfhnaf_dv(i) = dalphah_naf(k1) * (1.d0 - hnaf(i,L)) -
     X                  dbetah_naf(k1) * hnaf(i,L)
       dfmkdr_dv(i) = dalpham_kdr(k1) * (1.d0 - mkdr(i,L)) -
     X                  dbetam_kdr(k1) * mkdr(i,L)
       dfmka_dv(i)  = dalpham_ka(k1) * (1.d0 - mka(i,L)) -
     X                  dbetam_ka(k1) * mka(i,L)
       dfhka_dv(i)  = dalphah_ka(k1) * (1.d0 - hka(i,L)) -
     X                  dbetah_ka(k1) * hka(i,L)
       dfmk2_dv(i)  = dalpham_k2(k1) * (1.d0 - mk2(i,L)) -
     X                  dbetam_k2(k1) * mk2(i,L)
       dfhk2_dv(i)  = dalphah_k2(k1) * (1.d0 - hk2(i,L)) -
     X                  dbetah_k2(k1) * hk2(i,L)
       dfmkm_dv(i)  = dalpham_km(k1) * (1.d0 - mkm(i,L)) -
     X                  dbetam_km(k1) * mkm(i,L)
       dfmkc_dv(i)  = dalpham_kc(k1) * (1.d0 - mkc(i,L)) -
     X                  dbetam_kc(k1) * mkc(i,L)
       dfmcat_dv(i) = dalpham_cat(k1) * (1.d0 - mcat(i,L)) -
     X                  dbetam_cat(k1) * mcat(i,L)
       dfhcat_dv(i) = dalphah_cat(k1) * (1.d0 - hcat(i,L)) -
     X                  dbetah_cat(k1) * hcat(i,L)
       dfmcaL_dv(i) = dalpham_caL(k1) * (1.d0 - mcaL(i,L)) -
     X                  dbetam_caL(k1) * mcaL(i,L)
       dfmar_dv(i)  = dalpham_ar(k1) * (1.d0 - mar(i,L)) -
     X                  dbetam_ar(k1) * mar(i,L)

       dfmnaf_dmnaf(i) =  - alpham_naf(k0) - betam_naf(k0)
       dfmnap_dmnap(i) =  - alpham_naf(k2) - betam_naf(k2)
       dfhnaf_dhnaf(i) =  - alphah_naf(k1) - betah_naf(k1)
       dfmkdr_dmkdr(i) =  - alpham_kdr(k1) - betam_kdr(k1)
       dfmka_dmka(i)  =   - alpham_ka (k1) - betam_ka (k1)
       dfhka_dhka(i)  =   - alphah_ka (k1) - betah_ka (k1)
       dfmk2_dmk2(i)  =   - alpham_k2 (k1) - betam_k2 (k1)
       dfhk2_dhk2(i)  =   - alphah_k2 (k1) - betah_k2 (k1)
       dfmkm_dmkm(i)  =   - alpham_km (k1) - betam_km (k1)
       dfmkc_dmkc(i)  =   - alpham_kc (k1) - betam_kc (k1)
       dfmcat_dmcat(i) =  - alpham_cat(k1) - betam_cat(k1)
       dfhcat_dhcat(i) =  - alphah_cat(k1) - betah_cat(k1)
       dfmcaL_dmcaL(i) =  - alpham_caL(k1) - betam_caL(k1)
       dfmar_dmar(i)  =   - alpham_ar (k1) - betam_ar (k1)

          end do

       dt2 = 0.5d0 * dt * dt

        do i = 1, numcomp
          v(i,L) = v(i,L) + dt * fv(i)
           do j = 1, numcomp
        v(i,L) = v(i,L) + dt2 * dfv_dv(i,j) * fv(j)
           end do
        v(i,L) = v(i,L) + dt2 * ( dfv_dchi(i) * fchi(i)
     X          + dfv_dmnaf(i) * fmnaf(i)
     X          + dfv_dmnap(i) * fmnap(i)
     X          + dfv_dhnaf(i) * fhnaf(i)
     X          + dfv_dmkdr(i) * fmkdr(i)
     X          + dfv_dmka(i)  * fmka(i)
     X          + dfv_dhka(i)  * fhka(i)
     X          + dfv_dmk2(i)  * fmk2(i)
     X          + dfv_dhk2(i)  * fhk2(i)
     X          + dfv_dmkm(i)  * fmkm(i)
     X          + dfv_dmkc(i)  * fmkc(i)
     X          + dfv_dmkahp(i)* fmkahp(i)
     X          + dfv_dmcat(i)  * fmcat(i)
     X          + dfv_dhcat(i) * fhcat(i)
     X          + dfv_dmcaL(i) * fmcaL(i)
     X          + dfv_dmar(i)  * fmar(i) )

        chi(i,L) = chi(i,L) + dt * fchi(i) + dt2 *
     X   (dfchi_dchi(i) * fchi(i) + dfchi_dv(i) * fv(i))
        mnaf(i,L) = mnaf(i,L) + dt * fmnaf(i) + dt2 *
     X   (dfmnaf_dmnaf(i) * fmnaf(i) + dfmnaf_dv(i)*fv(i))
        mnap(i,L) = mnap(i,L) + dt * fmnap(i) + dt2 *
     X   (dfmnap_dmnap(i) * fmnap(i) + dfmnap_dv(i)*fv(i))
        hnaf(i,L) = hnaf(i,L) + dt * fhnaf(i) + dt2 *
     X   (dfhnaf_dhnaf(i) * fhnaf(i) + dfhnaf_dv(i)*fv(i))
        mkdr(i,L) = mkdr(i,L) + dt * fmkdr(i) + dt2 *
     X   (dfmkdr_dmkdr(i) * fmkdr(i) + dfmkdr_dv(i)*fv(i))
        mka(i,L) =  mka(i,L) + dt * fmka(i) + dt2 *
     X   (dfmka_dmka(i) * fmka(i) + dfmka_dv(i) * fv(i))
        hka(i,L) =  hka(i,L) + dt * fhka(i) + dt2 *
     X   (dfhka_dhka(i) * fhka(i) + dfhka_dv(i) * fv(i))
        mk2(i,L) =  mk2(i,L) + dt * fmk2(i) + dt2 *
     X   (dfmk2_dmk2(i) * fmk2(i) + dfmk2_dv(i) * fv(i))
        hk2(i,L) =  hk2(i,L) + dt * fhk2(i) + dt2 *
     X   (dfhk2_dhk2(i) * fhk2(i) + dfhk2_dv(i) * fv(i))
        mkm(i,L) =  mkm(i,L) + dt * fmkm(i) + dt2 *
     X   (dfmkm_dmkm(i) * fmkm(i) + dfmkm_dv(i) * fv(i))
        mkc(i,L) =  mkc(i,L) + dt * fmkc(i) + dt2 *
     X   (dfmkc_dmkc(i) * fmkc(i) + dfmkc_dv(i) * fv(i))
        mkahp(i,L) = mkahp(i,L) + dt * fmkahp(i) + dt2 *
     X (dfmkahp_dmkahp(i)*fmkahp(i) + dfmkahp_dchi(i)*fchi(i))
        mcat(i,L) =  mcat(i,L) + dt * fmcat(i) + dt2 *
     X   (dfmcat_dmcat(i) * fmcat(i) + dfmcat_dv(i) * fv(i))
        hcat(i,L) =  hcat(i,L) + dt * fhcat(i) + dt2 *
     X   (dfhcat_dhcat(i) * fhcat(i) + dfhcat_dv(i) * fv(i))
        mcaL(i,L) =  mcaL(i,L) + dt * fmcaL(i) + dt2 *
     X   (dfmcaL_dmcaL(i) * fmcaL(i) + dfmcaL_dv(i) * fv(i))
        mar(i,L) =   mar(i,L) + dt * fmar(i) + dt2 *
     X   (dfmar_dmar(i) * fmar(i) + dfmar_dv(i) * fv(i))
c            endif
         end do

! Add membrane currents into membcurr for appropriate compartments
          do i = 1, 9
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do
          do i = 14, 21
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do
          do i = 26, 33
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do
          do i = 39, 68
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do

            end do
c Finish loop L = 1 to numcell

         field_sup = 0.d0
         field_deep = 0.d0

         do i = 1, 12
        field_sup = field_sup + membcurr(i) / dabs(100.d0 - depth(i))
        field_deep = field_deep + membcurr(i) / dabs(500.d0 - depth(i))
         end do

2001          CONTINUE

6000    END



C  SETS UP TABLES FOR RATE FUNCTIONS
       SUBROUTINE SCORT_SETUP_L3pyr    
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)
      INTEGER I,J,K
      real*8 minf, hinf, taum, tauh, V, Z, shift_hnaf,
     X  shift_mkdr,
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
C FOR VOLTAGE, RANGE IS -120 TO +40 MV (absol.), 0.25 MV RESOLUTION


       DO 1, I = 0, 640
          V = dble(I)
          V = (V / 4.d0) - 120.d0

c gNa
           minf = 1.d0/(1.d0 + dexp((-V-38.d0)/10.d0))
           if (v.le.-30.d0) then
            taum = .025d0 + .14d0*dexp((v+30.d0)/10.d0)
           else
            taum = .02d0 + .145d0*dexp((-v-30.d0)/10.d0)
           endif
c from principal c. data, Martina & Jonas 1997, tau x 0.5
c Note that minf about the same for interneuron & princ. cell.
           alpham_naf(i) = minf / taum
           betam_naf(i) = 1.d0/taum - alpham_naf(i)

            shift_hnaf =  0.d0
        hinf = 1.d0/(1.d0 +
     x     dexp((v + shift_hnaf + 62.9d0)/10.7d0))
        tauh = 0.15d0 + 1.15d0/(1.d0+dexp((v+37.d0)/15.d0))
c from princ. cell data, Martina & Jonas 1997, tau x 0.5
            alphah_naf(i) = hinf / tauh
            betah_naf(i) = 1.d0/tauh - alphah_naf(i)

          shift_mkdr = 0.d0
c delayed rectifier, non-inactivating
       minf = 1.d0/(1.d0+dexp((-v-shift_mkdr-29.5d0)/10.0d0))
            if (v.le.-10.d0) then
             taum = .25d0 + 4.35d0*dexp((v+10.d0)/10.d0)
            else
             taum = .25d0 + 4.35d0*dexp((-v-10.d0)/10.d0)
            endif
              alpham_kdr(i) = minf / taum
              betam_kdr(i) = 1.d0 /taum - alpham_kdr(i)
c from Martina, Schultz et al., 1998. See espec. Table 1.

c A current: Huguenard & McCormick 1992, J Neurophysiol (TCR)
            minf = 1.d0/(1.d0 + dexp((-v-60.d0)/8.5d0))
            hinf = 1.d0/(1.d0 + dexp((v+78.d0)/6.d0))
        taum = .185d0 + .5d0/(dexp((v+35.8d0)/19.7d0) +
     x                            dexp((-v-79.7d0)/12.7d0))
        if (v.le.-63.d0) then
         tauh = .5d0/(dexp((v+46.d0)/5.d0) +
     x                  dexp((-v-238.d0)/37.5d0))
        else
         tauh = 9.5d0
        endif
           alpham_ka(i) = minf/taum
           betam_ka(i) = 1.d0 / taum - alpham_ka(i)
           alphah_ka(i) = hinf / tauh
           betah_ka(i) = 1.d0 / tauh - alphah_ka(i)

c h-current (anomalous rectifier), Huguenard & McCormick, 1992
           minf = 1.d0/(1.d0 + dexp((v+75.d0)/5.5d0))
           taum = 1.d0/(dexp(-14.6d0 -0.086d0*v) +
     x                   dexp(-1.87 + 0.07d0*v))
           alpham_ar(i) = minf / taum
           betam_ar(i) = 1.d0 / taum - alpham_ar(i)

c K2 K-current, McCormick & Huguenard
             minf = 1.d0/(1.d0 + dexp((-v-10.d0)/17.d0))
             hinf = 1.d0/(1.d0 + dexp((v+58.d0)/10.6d0))
            taum = 4.95d0 + 0.5d0/(dexp((v-81.d0)/25.6d0) +
     x                  dexp((-v-132.d0)/18.d0))
            tauh = 60.d0 + 0.5d0/(dexp((v-1.33d0)/200.d0) +
     x                  dexp((-v-130.d0)/7.1d0))
             alpham_k2(i) = minf / taum
             betam_k2(i) = 1.d0/taum - alpham_k2(i)
             alphah_k2(i) = hinf / tauh
             betah_k2(i) = 1.d0 / tauh - alphah_k2(i)

c voltage part of C-current, using 1994 kinetics, shift 60 mV
              if (v.le.-10.d0) then
       alpham_kc(i) = (2.d0/37.95d0)*dexp((v+50.d0)/11.d0 -
     x                                     (v+53.5)/27.d0)
       betam_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)-alpham_kc(i)
               else
       alpham_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)
       betam_kc(i) = 0.d0
               endif

c high-threshold gCa, from 1994, with 60 mV shift & no inactivn.
            alpham_cal(i) = 1.6d0/(1.d0+dexp(-.072d0*(v-5.d0)))
            betam_cal(i) = 0.1d0 * ((v+8.9d0)/5.d0) /
     x          (dexp((v+8.9d0)/5.d0) - 1.d0)

c M-current, from plast.f, with 60 mV shift
        alpham_km(i) = .02d0/(1.d0+dexp((-v-20.d0)/5.d0))
        betam_km(i) = .01d0 * dexp((-v-43.d0)/18.d0)

c T-current, from Destexhe, Neubig et al., 1998
         minf = 1.d0/(1.d0 + dexp((-v-56.d0)/6.2d0))
         hinf = 1.d0/(1.d0 + dexp((v+80.d0)/4.d0))
         taum = 0.204d0 + .333d0/(dexp((v+15.8d0)/18.2d0) +
     x                  dexp((-v-131.d0)/16.7d0))
          if (v.le.-81.d0) then
         tauh = 0.333 * dexp((v+466.d0)/66.6d0)
          else
         tauh = 9.32d0 + 0.333d0*dexp((-v-21.d0)/10.5d0)
          endif
              alpham_cat(i) = minf / taum
              betam_cat(i) = 1.d0/taum - alpham_cat(i)
              alphah_cat(i) = hinf / tauh
              betah_cat(i) = 1.d0 / tauh - alphah_cat(i)

1        CONTINUE

         do  i = 0, 639

      dalpham_naf(i) = (alpham_naf(i+1)-alpham_naf(i))/.25d0
      dbetam_naf(i) = (betam_naf(i+1)-betam_naf(i))/.25d0
      dalphah_naf(i) = (alphah_naf(i+1)-alphah_naf(i))/.25d0
      dbetah_naf(i) = (betah_naf(i+1)-betah_naf(i))/.25d0
      dalpham_kdr(i) = (alpham_kdr(i+1)-alpham_kdr(i))/.25d0
      dbetam_kdr(i) = (betam_kdr(i+1)-betam_kdr(i))/.25d0
      dalpham_ka(i) = (alpham_ka(i+1)-alpham_ka(i))/.25d0
      dbetam_ka(i) = (betam_ka(i+1)-betam_ka(i))/.25d0
      dalphah_ka(i) = (alphah_ka(i+1)-alphah_ka(i))/.25d0
      dbetah_ka(i) = (betah_ka(i+1)-betah_ka(i))/.25d0
      dalpham_k2(i) = (alpham_k2(i+1)-alpham_k2(i))/.25d0
      dbetam_k2(i) = (betam_k2(i+1)-betam_k2(i))/.25d0
      dalphah_k2(i) = (alphah_k2(i+1)-alphah_k2(i))/.25d0
      dbetah_k2(i) = (betah_k2(i+1)-betah_k2(i))/.25d0
      dalpham_km(i) = (alpham_km(i+1)-alpham_km(i))/.25d0
      dbetam_km(i) = (betam_km(i+1)-betam_km(i))/.25d0
      dalpham_kc(i) = (alpham_kc(i+1)-alpham_kc(i))/.25d0
      dbetam_kc(i) = (betam_kc(i+1)-betam_kc(i))/.25d0
      dalpham_cat(i) = (alpham_cat(i+1)-alpham_cat(i))/.25d0
      dbetam_cat(i) = (betam_cat(i+1)-betam_cat(i))/.25d0
      dalphah_cat(i) = (alphah_cat(i+1)-alphah_cat(i))/.25d0
      dbetah_cat(i) = (betah_cat(i+1)-betah_cat(i))/.25d0
      dalpham_caL(i) = (alpham_cal(i+1)-alpham_cal(i))/.25d0
      dbetam_caL(i) = (betam_cal(i+1)-betam_cal(i))/.25d0
      dalpham_ar(i) = (alpham_ar(i+1)-alpham_ar(i))/.25d0
      dbetam_ar(i) = (betam_ar(i+1)-betam_ar(i))/.25d0
       end do
2      CONTINUE

         do i = 640, 640
      dalpham_naf(i) =  dalpham_naf(i-1)
      dbetam_naf(i) =  dbetam_naf(i-1)
      dalphah_naf(i) = dalphah_naf(i-1)
      dbetah_naf(i) = dbetah_naf(i-1)
      dalpham_kdr(i) =  dalpham_kdr(i-1)
      dbetam_kdr(i) =  dbetam_kdr(i-1)
      dalpham_ka(i) =  dalpham_ka(i-1)
      dbetam_ka(i) =  dbetam_ka(i-1)
      dalphah_ka(i) =  dalphah_ka(i-1)
      dbetah_ka(i) =  dbetah_ka(i-1)
      dalpham_k2(i) =  dalpham_k2(i-1)
      dbetam_k2(i) =  dbetam_k2(i-1)
      dalphah_k2(i) =  dalphah_k2(i-1)
      dbetah_k2(i) =  dbetah_k2(i-1)
      dalpham_km(i) =  dalpham_km(i-1)
      dbetam_km(i) =  dbetam_km(i-1)
      dalpham_kc(i) =  dalpham_kc(i-1)
      dbetam_kc(i) =  dbetam_kc(i-1)
      dalpham_cat(i) =  dalpham_cat(i-1)
      dbetam_cat(i) =  dbetam_cat(i-1)
      dalphah_cat(i) =  dalphah_cat(i-1)
      dbetah_cat(i) =  dbetah_cat(i-1)
      dalpham_caL(i) =  dalpham_caL(i-1)
      dbetam_caL(i) =  dbetam_caL(i-1)
      dalpham_ar(i) =  dalpham_ar(i-1)
      dbetam_ar(i) =  dbetam_ar(i-1)
       end do   

4000   END

        SUBROUTINE SCORTMAJ_L3pyr   
C BRANCHED ACTIVE DENDRITES
     X             (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM,depth,level)
c Conductances: leak gL, coupling g, delayed rectifier gKDR, A gKA,
c C gKC, AHP gKAHP, K2 gK2, M gKM, low thresh Ca gCAT, high thresh
c gCAL, fast Na gNAF, persistent Na gNAP, h or anom. rectif. gAR.
c Note VAR = equil. potential for anomalous rectifier.
c Soma = comp. 1; 10 dendrites each with 13 compartments, 6-comp. axon
c Drop "glc"-like terms, just using "gl"-like
c CAFOR corresponds to "phi" in Traub et al., 1994
c Consistent set of units: nF, mV, ms, nA, microS

       INTEGER, PARAMETER:: numcomp = 74
! numcomp here must be compatible with numcomp_suppyrRS in calling prog.
        REAL*8 C(numcomp),GL(numcomp), GAM(0:numcomp, 0:numcomp)
        REAL*8 GNAF(numcomp),GCAT(numcomp), GKAHP(numcomp)
        REAL*8 GKDR(numcomp),GKA(numcomp),GKC(numcomp)
        REAL*8 GK2(numcomp),GNAP(numcomp),GAR(numcomp)
        REAL*8 GKM(numcomp), gcal(numcomp), CDENS
        REAL*8 JACOB(numcomp,numcomp),RI_SD,RI_AXON,RM_SD,RM_AXON
        INTEGER LEVEL(numcomp)
        REAL*8 GNAF_DENS(0:12), GCAT_DENS(0:12), GKDR_DENS(0:12)
        REAL*8 GKA_DENS(0:12), GKC_DENS(0:12), GKAHP_DENS(0:12)
        REAL*8 GCAL_DENS(0:12), GK2_DENS(0:12), GKM_DENS(0:12)
        REAL*8 GNAP_DENS(0:12), GAR_DENS(0:12)
        REAL*8 RES, RINPUT, Z, ELEN(numcomp)
        REAL*8 RSOMA, PI, BETCHI(numcomp), CAFOR(numcomp)
        REAL*8 RAD(numcomp), LEN(numcomp), GAM1, GAM2
        REAL*8 RIN, D(numcomp), AREA(numcomp), RI
        INTEGER NEIGH(numcomp,10), NNUM(numcomp), i, j, k, it
C FOR ESTABLISHING TOPOLOGY OF COMPARTMENTS
        real*8 depth(12) ! depth in microns of levels 1-12, assuming soma 
! at depth 500 microns 

        depth(1) = 600.d0
        depth(2) = 650.d0
        depth(3) = 700.d0
        depth(4) = 750.d0
        depth(5) = 550.d0
        depth(6) = 450.d0
        depth(7) = 400.d0
        depth(8) = 350.d0
        depth(9) = 300.d0
        depth(10) = 250.d0
        depth(11) = 200.d0
        depth(12) =  50.d0

        RI_SD = 250.d0
        RM_SD = 50000.d0
        RI_AXON = 100.d0
        RM_AXON = 1000.d0
        CDENS = 0.9d0

        PI = 3.14159d0

       do i = 0, 12
        gnaf_dens(i) = 10.d0
       end do
c       gnaf_dens(0) = 400.d0
!       gnaf_dens(0) = 120.d0
        gnaf_dens(0) = 200.d0
        gnaf_dens(1) = 120.d0
        gnaf_dens(2) =  75.d0
        gnaf_dens(5) = 100.d0
        gnaf_dens(6) =  75.d0

       do i = 0, 12
        gkdr_dens(i) = 0.d0
       end do
c       gkdr_dens(0) = 400.d0
c       gkdr_dens(0) = 100.d0
c       gkdr_dens(0) = 170.d0
        gkdr_dens(0) = 250.d0
c       gkdr_dens(1) = 100.d0
        gkdr_dens(1) = 150.d0
        gkdr_dens(2) =  75.d0
        gkdr_dens(5) = 100.d0
        gkdr_dens(6) =  75.d0

        gnap_dens(0) = 0.d0
        do i = 1, 12
          gnap_dens(i) = 0.0040d0 * gnaf_dens(i)
c         gnap_dens(i) = 0.002d0 * gnaf_dens(i)
c         gnap_dens(i) = 0.0030d0 * gnaf_dens(i)
        end do

        gcat_dens(0) = 0.d0
        do i = 1, 12
c         gcat_dens(i) = 0.5d0
          gcat_dens(i) = 0.1d0
        end do

        gcaL_dens(0) = 0.d0
        do i = 1, 6
          gcaL_dens(i) = 0.5d0
        end do
        do i = 7, 12
          gcaL_dens(i) = 0.5d0
        end do

       do i = 0, 12
        gka_dens(i) = 2.d0
       end do
        gka_dens(0) =100.d0 ! NOTE
        gka_dens(1) = 30.d0
        gka_dens(5) = 30.d0

      do i = 0, 12
c        gkc_dens(i)  = 12.00d0
         gkc_dens(i)  =  0.00d0
c        gkc_dens(i)  =  2.00d0
c        gkc_dens(i)  =  7.00d0
      end do
         gkc_dens(0) =  0.00d0
c        gkc_dens(1) = 7.5d0
c        gkc_dens(1) = 12.d0
         gkc_dens(1) = 15.d0
c        gkc_dens(2) = 7.5d0
         gkc_dens(2) = 10.d0
         gkc_dens(5) = 7.5d0
         gkc_dens(6) = 7.5d0

c       gkm_dens(0) = 2.d0 ! 9 Nov. 2005, see scort-pan.f of today
        gkm_dens(0) = 8.d0 ! 9 Nov. 2005, see scort-pan.f of today
! Above suppresses doublets, but still allows FRB with appropriate
! gNaP, gKC, and rel_axonshift (e.g. 6 mV)
        do i = 1, 12
         gkm_dens(i) = 2.5d0 * 1.50d0
        end do

        do i = 0, 12
c       gk2_dens(i) = 1.d0
        gk2_dens(i) = 0.1d0
        end do
        gk2_dens(0) = 0.d0

        gkahp_dens(0) = 0.d0
        do i = 1, 12
c        gkahp_dens(i) = 0.200d0
         gkahp_dens(i) = 0.100d0
c        gkahp_dens(i) = 0.050d0
        end do

        gar_dens(0) = 0.d0
        do i = 1, 12
         gar_dens(i) = 0.25d0
        end do

c       WRITE   (6,9988)
9988    FORMAT(2X,'I',4X,'NADENS',' CADENS(T)',' KDRDEN',' KAHPDE',
     X     ' KCDENS',' KADENS')
        DO 9989, I = 0, 12
c         WRITE (6,9990) I, gnaf_dens(i), gcat_dens(i), gkdr_dens(i),
c    X  gkahp_dens(i), gkc_dens(i), gka_dens(i)
9990    FORMAT(2X,I2,2X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2)
9989    CONTINUE


        level(1) = 1
        do i = 2, 13
         level(i) = 2
        end do
        do i = 14, 25
           level(i) = 3
        end do
        do i = 26, 37
           level(i) = 4
        end do
        level(38) = 5
        level(39) = 6
        level(40) = 7
        level(41) = 8
        level(42) = 8
        level(43) = 9
        level(44) = 9
        do i = 45, 52
           level(i) = 10
        end do
        do i = 53, 60
           level(i) = 11
        end do
        do i = 61, 68
           level(i) = 12
        end do

        do i =  69, 74
         level(i) = 0
        end do

c connectivity of axon
        nnum( 69) = 2
        nnum( 70) = 3
        nnum( 71) = 3
        nnum( 73) = 3
        nnum( 72) = 1
        nnum( 74) = 1
         neigh(69,1) =  1
         neigh(69,2) = 70
         neigh(70,1) = 69
         neigh(70,2) = 71
         neigh(70,3) = 73
         neigh(71,1) = 70
         neigh(71,2) = 72
         neigh(71,3) = 73
         neigh(73,1) = 70
         neigh(73,2) = 71
         neigh(73,3) = 74
         neigh(72,1) = 71
         neigh(74,1) = 73

c connectivity of SD part
          nnum(1) = 10
          neigh(1,1) = 69
          neigh(1,2) =  2
          neigh(1,3) =  3
          neigh(1,4) =  4
          neigh(1,5) =  5
          neigh(1,6) =  6
          neigh(1,7) =  7
          neigh(1,8) =  8
          neigh(1,9) =  9
          neigh(1,10) = 38

          do i = 2, 9
           nnum(i) = 2
           neigh(i,1) = 1
           neigh(i,2) = i + 12
          end do

          do i = 14, 21
            nnum(i) = 2
            neigh(i,1) = i - 12
            neigh(i,2) = i + 12
          end do

          do i = 26, 33
            nnum(i) = 1
            neigh(i,1) = i - 12
          end do

          do i = 10, 13
            nnum(i) = 2
            neigh(i,1) = 38
            neigh(i,2) = i + 12
          end do

          do i = 22, 25
            nnum(i) = 2
            neigh(i,1) = i - 12
            neigh(i,2) = i + 12
          end do

          do i = 34, 37
            nnum(i) = 1
            neigh(i,1) = i - 12
          end do

          nnum(38) = 6
          neigh(38,1) = 1
          neigh(38,2) = 39
          neigh(38,3) = 10
          neigh(38,4) = 11
          neigh(38,5) = 12
          neigh(38,6) = 13

          nnum(39) = 2
          neigh(39,1) = 38
          neigh(39,2) = 40

          nnum(40) = 3
          neigh(40,1) = 39
          neigh(40,2) = 41
          neigh(40,3) = 42

          nnum(41) = 3
          neigh(41,1) = 40
          neigh(41,2) = 42
          neigh(41,3) = 43

          nnum(42) = 3
          neigh(42,1) = 40
          neigh(42,2) = 41
          neigh(42,3) = 44

           nnum(43) = 5
           neigh(43,1) = 41
           neigh(43,2) = 45
           neigh(43,3) = 46
           neigh(43,4) = 47
           neigh(43,5) = 48

           nnum(44) = 5
           neigh(44,1) = 42
           neigh(44,2) = 49
           neigh(44,3) = 50
           neigh(44,4) = 51
           neigh(44,5) = 52

           nnum(45) = 5
           neigh(45,1) = 43
           neigh(45,2) = 53
           neigh(45,3) = 46
           neigh(45,4) = 47
           neigh(45,5) = 48

           nnum(46) = 5
           neigh(46,1) = 43
           neigh(46,2) = 54
           neigh(46,3) = 45
           neigh(46,4) = 47
           neigh(46,5) = 48

           nnum(47) = 5
           neigh(47,1) = 43
           neigh(47,2) = 55
           neigh(47,3) = 45
           neigh(47,4) = 46
           neigh(47,5) = 48

           nnum(48) = 5
           neigh(48,1) = 43
           neigh(48,2) = 56
           neigh(48,3) = 45
           neigh(48,4) = 46
           neigh(48,5) = 47

           nnum(49) = 5
           neigh(49,1) = 44
           neigh(49,2) = 57
           neigh(49,3) = 50
           neigh(49,4) = 51
           neigh(49,5) = 52

           nnum(50) = 5
           neigh(50,1) = 44
           neigh(50,2) = 58
           neigh(50,3) = 49
           neigh(50,4) = 51
           neigh(50,5) = 52

           nnum(51) = 5
           neigh(51,1) = 44
           neigh(51,2) = 59
           neigh(51,3) = 49
           neigh(51,4) = 50
           neigh(51,5) = 52

           nnum(52) = 5
           neigh(52,1) = 44
           neigh(52,2) = 60
           neigh(52,3) = 49
           neigh(52,4) = 51
           neigh(52,5) = 50

          do i = 53, 60
           nnum(i) = 2
           neigh(i,1) = i - 8
           neigh(i,2) = i + 8
          end do

          do i = 61, 68
           nnum(i) = 1
           neigh(i,1) = i - 8
          end do

c        DO 332, I = 1, 74
         DO I = 1, 74
c          WRITE(6,3330) I, NEIGH(I,1),NEIGH(I,2),NEIGH(I,3),NEIGH(I,4),
c    X NEIGH(I,5),NEIGH(I,6),NEIGH(I,7),NEIGH(I,8),NEIGH(I,9),
c    X NEIGH(I,10)
3330     FORMAT(2X,11I5)
         END DO
332      CONTINUE
c         DO 858, I = 1, 74
          DO I = 1, 74
c          DO 858, J = 1, NNUM(I)
           DO J = 1, NNUM(I)
            K = NEIGH(I,J)
            IT = 0
c           DO 859, L = 1, NNUM(K)
            DO  L = 1, NNUM(K)
             IF (NEIGH(K,L).EQ.I) IT = 1
            END DO
859         CONTINUE
             IF (IT.EQ.0) THEN
c             WRITE(6,8591) I, K
8591          FORMAT(' ASYMMETRY IN NEIGH MATRIX ',I4,I4)
              STOP
             ENDIF
          END DO
          END DO
858       CONTINUE

c length and radius of axonal compartments
c Note shortened "initial segment"
          len(69) = 25.d0
          do i = 70, 74
            len(i) = 50.d0
          end do
          rad( 69) = 0.90d0
c         rad( 69) = 0.80d0
          rad( 70) = 0.7d0
          do i = 71, 74
           rad(i) = 0.5d0
          end do

c  length and radius of SD compartments
          len(1) = 15.d0
          rad(1) =  8.d0

          do i = 2, 68
           len(i) = 50.d0
          end do
c lengthen some compartments, e.g. apical shaft
           len(38) = 65.d0
           len(39) = 65.d0
           len(40) = 65.d0

          do i = 2, 37
            rad(i) = 0.5d0
          end do

          z = 4.0d0
          rad(38) = z
          rad(39) = 0.9d0 * z
          rad(40) = 0.8d0 * z
          rad(41) = 0.5d0 * z
          rad(42) = 0.5d0 * z
          rad(43) = 0.5d0 * z
          rad(44) = 0.5d0 * z
          do i = 45, 68
           rad(i) = 0.2d0 * z
          end do


c       WRITE(6,919)
919     FORMAT('COMPART.',' LEVEL ',' RADIUS ',' LENGTH(MU)')
c       DO 920, I = 1, 74
c920      WRITE(6,921) I, LEVEL(I), RAD(I), LEN(I)
921     FORMAT(I3,5X,I2,3X,F6.2,1X,F6.1,2X,F4.3)

        DO 120, I = 1, 74
          AREA(I) = 2.d0 * PI * RAD(I) * LEN(I)
      if((i.gt.1).and.(i.le.68)) area(i) = 2.d0 * area(i)
C    CORRECTION FOR CONTRIBUTION OF SPINES TO AREA
          K = LEVEL(I)
          C(I) = CDENS * AREA(I) * (1.D-8)

           if (k.ge.1) then
          GL(I) = (1.D-2) * AREA(I) / RM_SD
           else
          GL(I) = (1.D-2) * AREA(I) / RM_AXON
           endif

          GNAF(I) = GNAF_DENS(K) * AREA(I) * (1.D-5)
          GNAP(I) = GNAP_DENS(K) * AREA(I) * (1.D-5)
          GCAT(I) = GCAT_DENS(K) * AREA(I) * (1.D-5)
          GKDR(I) = GKDR_DENS(K) * AREA(I) * (1.D-5)
          GKA(I) = GKA_DENS(K) * AREA(I) * (1.D-5)
          GKC(I) = GKC_DENS(K) * AREA(I) * (1.D-5)
          GKAHP(I) = GKAHP_DENS(K) * AREA(I) * (1.D-5)
          GKM(I) = GKM_DENS(K) * AREA(I) * (1.D-5)
          GCAL(I) = GCAL_DENS(K) * AREA(I) * (1.D-5)
          GK2(I) = GK2_DENS(K) * AREA(I) * (1.D-5)
          GAR(I) = GAR_DENS(K) * AREA(I) * (1.D-5)
c above conductances should be in microS
120           continue

         Z = 0.d0
c        DO 1019, I = 2, 68
         DO I = 2, 68
           Z = Z + AREA(I)
         END DO
1019     CONTINUE
c        WRITE(6,1020) Z
1020     FORMAT(2X,' TOTAL DENDRITIC AREA ',F7.0)

c       DO 140, I = 1, 74
        DO I = 1, 74
c       DO 140, K = 1, NNUM(I)
        DO K = 1, NNUM(I)
         J = NEIGH(I,K)
           if (level(i).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM1 =100.d0 * PI * RAD(I) * RAD(I) / ( RI * LEN(I) )

           if (level(j).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM2 =100.d0 * PI * RAD(J) * RAD(J) / ( RI * LEN(J) )
         GAM(I,J) = 2.d0/( (1.d0/GAM1) + (1.d0/GAM2) )
	 END DO
	 END DO

140     CONTINUE
c gam computed in microS

c       DO 299, I = 1, 74
        DO I = 1, 74
299       BETCHI(I) = .05d0
        END DO
        BETCHI( 1) =  .01d0

c       DO 300, I = 1, 74
        DO I = 1, 74
c300     D(I) = 2.D-4
300     D(I) = 5.D-4
        END DO
c       DO 301, I = 1, 74
        DO I = 1, 74
         IF (LEVEL(I).EQ.1) D(I) = 2.D-3
        END DO
301     CONTINUE
C  NOTE NOTE NOTE  (DIFFERENT FROM SWONG)


c      DO 160, I = 1, 74
       DO I = 1, 74
160     CAFOR(I) = 5200.d0 / (AREA(I) * D(I))
       END DO
C     NOTE CORRECTION

c       do 200, i = 1, 74
        do i = 1, numcomp
200     C(I) = 1000.d0 * C(I)
        end do
C     TO GO FROM MICROF TO NF.

c     DO 909, I = 1, 74
      DO I = 1, numcomp
       JACOB(I,I) = - GL(I)
c     DO 909, J = 1, NNUM(I)
      DO J = 1, NNUM(I)
         K = NEIGH(I,J)
         IF (I.EQ.K) THEN
c            WRITE(6,510) I
510          FORMAT(' UNEXPECTED SYMMETRY IN NEIGH ',I4)
         ENDIF
         JACOB(I,K) = GAM(I,K)
         JACOB(I,I) = JACOB(I,I) - GAM(I,K)
       END DO
       END DO
909   CONTINUE

c 15 Jan. 2001: make correction for c(i)
          do i = 1, numcomp
          do j = 1, numcomp
             jacob(i,j) = jacob(i,j) / c(i)
          end do
          end do

c      DO 500, I = 1, 74
       DO I = 1, 74
c       WRITE (6,501) I,C(I)
501     FORMAT(1X,I3,' C(I) = ',F7.4)
       END DO
500     CONTINUE
        END


c 22 Aug 2019, start with suppyrRS integration subroutine from
c son_of_groucho, and use for semilunar in piriform simulations.
c Need to change field variables and depth definitions,
c  and perhaps alter compartment dimensions.
c Also disconnect basal dendrites.

c 11 Sept 2006, start with /interact/integrate_suppyrRSXP.f & add GABA-B
! 7 Nov. 2005: modify integrate_suppyrRSX.f to allow for Colbert-Pan axon.
!29 July 2005: modify groucho/integrate_suppyrRS.f, for a separate
! call for initialization, and to integrate only selected cells.
! Integration routine for suppyrRS cells
! Routine adapted from scortn in supergj.f
c      SUBROUTINE INTEGRATE_suppyrRSXPB (O, time, numcell,     
       SUBROUTINE INTEGRATE_semilunar   (O, time, numcell,     
     &    V, curr, initialize, firstcell, lastcell,
     & gAMPA, gNMDA, gGABA_A, gGABA_B,
     & Mg, 
     & gapcon  ,totaxgj   ,gjtable, dt,
     &  chi,mnaf,mnap,
     &  hnaf,mkdr,mka,
     &  hka,mk2,hk2,
     &  mkm,mkc,mkahp,
     &  mcat,hcat,mcal,
     &  mar,field_sup,field_deep,rel_axonshift)

       SAVE

       INTEGER, PARAMETER:: numcomp = 74
! numcomp here must be compatible with numcomp_suppyrRS in calling prog.
       INTEGER  numcell, num_other
       INTEGER initialize, firstcell, lastcell
       INTEGER J1, I, J, K, K1, K2, K3, L, L1, O
       REAL*8 c(numcomp), curr(numcomp,numcell)
       REAL*8  Z, Z1, Z2, Z3, Z4, DT, time
       integer totaxgj, gjtable(totaxgj,4)
       real*8 gapcon, gAMPA(numcomp,numcell),
     &        gNMDA(numcomp,numcell), gGABA_A(numcomp,numcell),
     &        gGABA_B(numcomp,numcell)
       real*8 Mg, V(numcomp,numcell), rel_axonshift

c CINV is 1/C, i.e. inverse capacitance
       real*8 chi(numcomp,numcell),
     & mnaf(numcomp,numcell),mnap(numcomp,numcell),
     x hnaf(numcomp,numcell), mkdr(numcomp,numcell),
     x mka(numcomp,numcell),hka(numcomp,numcell),
     x mk2(numcomp,numcell), cinv(numcomp),
     x hk2(numcomp,numcell),mkm(numcomp,numcell),
     x mkc(numcomp,numcell),mkahp(numcomp,numcell),
     x mcat(numcomp,numcell),hcat(numcomp,numcell),
     x mcal(numcomp,numcell), betchi(numcomp),
     x mar(numcomp,numcell),jacob(numcomp,numcomp),
     x gam(0: numcomp,0: numcomp),gL(numcomp),gnaf(numcomp),
     x gnap(numcomp),gkdr(numcomp),gka(numcomp),
     x gk2(numcomp),gkm(numcomp),
     x gkc(numcomp),gkahp(numcomp),
     x gcat(numcomp),gcaL(numcomp),gar(numcomp),
     x cafor(numcomp)
       real*8
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
       real*8 vL(numcomp),vk(numcomp),vna,var,vca,vgaba_a
       real*8 depth(12), membcurr(12), field_sup, field_deep
       integer level(numcomp)

        INTEGER NEIGH(numcomp,10), NNUM(numcomp)
        INTEGER igap1, igap2
c the f's are the functions giving 1st derivatives for evolution of
c the differential equations for the voltages (v), calcium (chi), and
c other state variables.
       real*8 fv(numcomp), fchi(numcomp),
     x fmnaf(numcomp),fhnaf(numcomp),fmkdr(numcomp),
     x fmka(numcomp),fhka(numcomp),fmk2(numcomp),
     x fhk2(numcomp),fmnap(numcomp),
     x fmkm(numcomp),fmkc(numcomp),fmkahp(numcomp),
     x fmcat(numcomp),fhcat(numcomp),fmcal(numcomp),
     x fmar(numcomp)

c below are for calculating the partial derivatives
       real*8 dfv_dv(numcomp,numcomp), dfv_dchi(numcomp),
     x  dfv_dmnaf(numcomp),  dfv_dmnap(numcomp),
     x  dfv_dhnaf(numcomp),dfv_dmkdr(numcomp),
     x  dfv_dmka(numcomp),dfv_dhka(numcomp),
     x  dfv_dmk2(numcomp),dfv_dhk2(numcomp),
     x  dfv_dmkm(numcomp),dfv_dmkc(numcomp),
     x  dfv_dmkahp(numcomp),dfv_dmcat(numcomp),
     x  dfv_dhcat(numcomp),dfv_dmcal(numcomp),
     x  dfv_dmar(numcomp)

        real*8 dfchi_dv(numcomp), dfchi_dchi(numcomp),
     x dfmnaf_dmnaf(numcomp), dfmnaf_dv(numcomp),
     x dfhnaf_dhnaf(numcomp),
     x dfmnap_dmnap(numcomp), dfmnap_dv(numcomp),
     x dfhnaf_dv(numcomp),dfmkdr_dmkdr(numcomp),
     x dfmkdr_dv(numcomp),
     x dfmka_dmka(numcomp),dfmka_dv(numcomp),
     x dfhka_dhka(numcomp),dfhka_dv(numcomp),
     x dfmk2_dmk2(numcomp),dfmk2_dv(numcomp),
     x dfhk2_dhk2(numcomp),dfhk2_dv(numcomp),
     x dfmkm_dmkm(numcomp),dfmkm_dv(numcomp),
     x dfmkc_dmkc(numcomp),dfmkc_dv(numcomp),
     x dfmcat_dmcat(numcomp),dfmcat_dv(numcomp),dfhcat_dhcat(numcomp),
     x dfhcat_dv(numcomp),dfmcal_dmcal(numcomp),dfmcal_dv(numcomp),
     x dfmar_dmar(numcomp),dfmar_dv(numcomp),dfmkahp_dchi(numcomp),
     x dfmkahp_dmkahp(numcomp), dt2

       REAL*8 OPEN(numcomp),gamma(numcomp),gamma_prime(numcomp)
c gamma is function of chi used in calculating KC conductance
       REAL*8 alpham_ahp(numcomp), alpham_ahp_prime(numcomp)
       REAL*8 gna_tot(numcomp),gk_tot(numcomp),gca_tot(numcomp)
       REAL*8 gca_high(numcomp), gar_tot(numcomp)
c this will be gCa conductance corresponding to high-thresh channels

       real*8 persistentNa_shift, fastNa_shift_SD,
     x   fastNa_shift_axon

       REAL*8 A, BB1, BB2  ! params. for FNMDA.f


c          if (O.eq.1) then
           if (initialize.eq.0) then
c do initialization

c Program fnmda assumes A, BB1, BB2 defined in calling program
c as follows:
         A = DEXP(-2.847d0)
         BB1 = DEXP(-.693d0)
         BB2 = DEXP(-3.101d0)

c       goto 4000
       CALL   SCORT_SETUP_semilunar
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)

        CALL SCORTMAJ_semilunar
     X             (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM,depth,level)

          do i = 1, numcomp
             cinv(i) = 1.d0 / c(i)
          end do
4000      CONTINUE

           do i = 1, numcomp
          vL(i) = -70.d0
          vK(i) = -95.d0
           end do

        VNA = 50.d0
        VCA = 125.d0
        VAR = -43.d0
        VAR = -35.d0
c -43 mV from Huguenard & McCormick
        VGABA_A = -81.d0
c       write(6,901) VNa, VCa, VK(1), O
901     format('VNa =',f6.2,' VCa =',f6.2,' VK =',f6.2,
     &   ' O = ',i3)

c ? initialize membrane state variables?
         do L = 1, numcell  
         do i = 1, numcomp
        v(i,L) = VL(i)
	chi(i,L) = 0.d0
	mnaf(i,L) = 0.d0
	mkdr(i,L) = 0.d0
	mk2(i,L) = 0.d0
	mkm(i,L) = 0.d0
	mkc(i,L) = 0.d0
	mkahp(i,L) = 0.d0
	mcat(i,L) = 0.d0
	mcal(i,L) = 0.d0
         end do
         end do

          do L = 1, numcell
        k1 = idnint (4.d0 * (v(1,L) + 120.d0))

            do i = 1, numcomp
      hnaf(i,L) = alphah_naf(k1)/(alphah_naf(k1)
     &       +betah_naf(k1))
      hka(i,L) = alphah_ka(k1)/(alphah_ka(k1)
     &                               +betah_ka(k1))
      hk2(i,L) = alphah_k2(k1)/(alphah_k2(k1)
     &                                +betah_k2(k1))
      hcat(i,L)=alphah_cat(k1)/(alphah_cat(k1)
     &                                +betah_cat(k1))
c     mar=alpham_ar(k1)/(alpham_ar(k1)+betam_ar(k1))
      mar(i,L) = .25d0
             end do
           end do


             do i = 1, numcomp
	    open(i) = 0.d0
            gkm(i) = 2.d0 * gkm(i)
             end do

         do i = 1, 68
c          gnaf(i) = 0.8d0 * 1.25d0 * gnaf(i) ! factor of 0.8 added 19 Nov. 2005
c          gnaf(i) = 0.9d0 * 1.25d0 * gnaf(i) ! Back to 0.9, 29 Nov. 2005
           gnaf(i) = 0.6d0 * 1.25d0 * gnaf(i) ! 
! NOTE THAT THERE IS QUESTION OF HOW TO COMPARE BEHAVIOR OF PYRAMID IN NETWORK WITH
! SIMULATIONS OF SINGLE CELL.  IN FORMER CASE, THERE IS LARGE AXONAL SHUNT THROUGH
! gj(s), NOT PRESENT IN SINGLE CELL MODEL.  THEREFORE, HIGHER AXONAL gNa MIGHT BE
! NECESSARY FOR SPIKE PROPAGATION.
c          gnaf(i) = 0.9d0 * 1.25d0 * gnaf(i) ! factor of 0.9 added 20 Nov. 2005
           gkdr(i) = 1.25d0 * gkdr(i)
         end do
 
c Perhaps reduce fast gNa on IS
          gnaf(69) = 1.00d0 * gnaf(69)
c         gnaf(69) = 0.25d0 * gnaf(69)
          gnaf(70) = 1.00d0 * gnaf(70)
c         gnaf(70) = 0.25d0 * gnaf(70)

c Perhaps reduce coupling between soma and IS
c         gam(1,69) = 0.15d0 * gam(1,69)
c         gam(69,1) = 0.15d0 * gam(69,1)

               z1 = 0.0d0
c              z2 = 1.2d0 ! value 1.2 tried Feb. 21, 2013
               z2 = 1.5d0 ! value 1.2 tried Feb. 21, 2013
               z3 = 1.0d0
c              z3 = 0.0d0 ! Note reduction from 0.4, to prevent
c slow hyperpolarization that seems to mess up gamma.
               z4 = 0.3d0
c RS cell
             do i = 1, numcomp
              gnap(i) = z1 * gnap(i)
              gkc (i) = z2 * gkc (i)
              gkahp(i) = z3 * gkahp(i)
              gkm (i) = z4 * gkm (i)
             end do

              goto 6000

          endif
c End initialization

          do i = 1, 12
           membcurr(i) = 0.d0
          end do

c                  goto 2001


c             do L = 1, numcell
              do L = firstcell, lastcell

	  do i = 1, numcomp
	  do j = 1, nnum(i)
	   if (neigh(i,j).gt.numcomp) then
          write(6,433) i, j, L
433       format(' ls ',3x,3i5)
           endif
	end do
	end do

       DO I = 1, numcomp
          FV(I) = -GL(I) * (V(I,L) - VL(i)) * cinv(i)
          DO J = 1, NNUM(I)
             K = NEIGH(I,J)
302     FV(I) = FV(I) + GAM(I,K) * (V(K,L) - V(I,L)) * cinv(i)
           END DO
       END DO
301    CONTINUE


       CALL FNMDA (V, OPEN, numcell, numcomp, MG, L,
     &                 A, BB1, BB2)

      DO I = 1, numcomp
       FV(I) = FV(I) + ( CURR(I,L)
     X   - (gampa(I,L) + open(i) * gnmda(I,L))*V(I,L)
     X   - ggaba_a(I,L)*(V(I,L)-Vgaba_a) 
     X   - ggaba_b(I,L)*(V(I,L)-VK(i)  ) ) * cinv(i)
c above assumes equil. potential for AMPA & NMDA = 0 mV
      END DO
421      continue

       do m = 1, totaxgj
        if (gjtable(m,1).eq.L) then
         L1 = gjtable(m,3)
         igap1 = gjtable(m,2)
         igap2 = gjtable(m,4)
 	fv(igap1) = fv(igap1) + gapcon *
     &   (v(igap2,L1) - v(igap1,L)) * cinv(igap1)
        else if (gjtable(m,3).eq.L) then
         L1 = gjtable(m,1)
         igap1 = gjtable(m,4)
         igap2 = gjtable(m,2)
 	fv(igap1) = fv(igap1) + gapcon *
     &   (v(igap2,L1) - v(igap1,L)) * cinv(igap1)
        endif
       end do ! do m


       do i = 1, numcomp
        gamma(i) = dmin1 (1.d0, .004d0 * chi(i,L))
        if (chi(i,L).le.250.d0) then
          gamma_prime(i) = .004d0
        else
          gamma_prime(i) = 0.d0
        endif
c         endif
       end do

      DO I = 1, numcomp
       gna_tot(i) = gnaf(i) * (mnaf(i,L)**3) * hnaf(i,L) +
     x     gnap(i) * mnap(i,L)
       gk_tot(i) = gkdr(i) * (mkdr(i,L)**4) +
     x             gka(i)  * (mka(i,L)**4) * hka(i,L) +
     x             gk2(i)  * mk2(i,L) * hk2(i,L) +
     x             gkm(i)  * mkm(i,L) +
     x             gkc(i)  * mkc(i,L) * gamma(i) +
     x             gkahp(i)* mkahp(i,L)
       gca_tot(i) = gcat(i) * (mcat(i,L)**2) * hcat(i,L) +
     x              gcaL(i) * (mcaL(i,L)**2)
       gca_high(i) =
     x              gcaL(i) * (mcaL(i,L)**2)
       gar_tot(i) = gar(i) * mar(i,L)


       FV(I) = FV(I) - ( gna_tot(i) * (v(i,L) - vna)
     X  + gk_tot(i) * (v(i,L) - vK(i))
     X  + gca_tot(i) * (v(i,L) - vCa)
     X  + gar_tot(i) * (v(i,L) - var) ) * cinv(i)
c        endif
       END DO
88           continue

         do i = 1, numcomp
         do j = 1, numcomp
          if (i.ne.j) then
            dfv_dv(i,j) = jacob(i,j)
          else
            dfv_dv(i,j) = jacob(i,i) - cinv(i) *
     X  (gna_tot(i) + gk_tot(i) + gca_tot(i) + gar_tot(i)
     X   + ggaba_a(i,L) + ggaba_b(i,L) + gampa(i,L)
     X   + open(i) * gnmda(I,L) )
          endif
         end do
         end do

           do i = 1, numcomp
        dfv_dchi(i)  = - cinv(i) * gkc(i) * mkc(i,L) *
     x                     gamma_prime(i) * (v(i,L)-vK(i))
        dfv_dmnaf(i) = -3.d0 * cinv(i) * (mnaf(i,L)**2) *
     X    (gnaf(i) * hnaf(i,L)          ) * (v(i,L) - vna)
        dfv_dmnap(i) = - cinv(i) *
     X    (               gnap(i)) * (v(i,L) - vna)
        dfv_dhnaf(i) = - cinv(i) * gnaf(i) * (mnaf(i,L)**3) *
     X                    (v(i,L) - vna)
        dfv_dmkdr(i) = -4.d0 * cinv(i) * gkdr(i) * (mkdr(i,L)**3)
     X                   * (v(i,L) - vK(i))
        dfv_dmka(i)  = -4.d0 * cinv(i) * gka(i) * (mka(i,L)**3) *
     X                   hka(i,L) * (v(i,L) - vK(i))
        dfv_dhka(i)  = - cinv(i) * gka(i) * (mka(i,L)**4) *
     X                    (v(i,L) - vK(i))
      dfv_dmk2(i) = - cinv(i) * gk2(i) * hk2(i,L) * (v(i,L)-vK(i))
      dfv_dhk2(i) = - cinv(i) * gk2(i) * mk2(i,L) * (v(i,L)-vK(i))
      dfv_dmkm(i) = - cinv(i) * gkm(i) * (v(i,L) - vK(i))
      dfv_dmkc(i) = - cinv(i)*gkc(i) * gamma(i) * (v(i,L)-vK(i))
        dfv_dmkahp(i)= - cinv(i) * gkahp(i) * (v(i,L) - vK(i))
        dfv_dmcat(i)  = -2.d0 * cinv(i) * gcat(i) * mcat(i,L) *
     X                    hcat(i,L) * (v(i,L) - vCa)
        dfv_dhcat(i) = - cinv(i) * gcat(i) * (mcat(i,L)**2) *
     X                  (v(i,L) - vCa)
        dfv_dmcal(i) = -2.d0 * cinv(i) * gcal(i) * mcal(i,L) *
     X                      (v(i,L) - vCa)
        dfv_dmar(i) = - cinv(i) * gar(i) * (v(i,L) - var)
            end do

         do i = 1, numcomp
          fchi(i) = - cafor(i) * gca_high(i) * (v(i,L) - vca)
     x       - betchi(i) * chi(i,L)
          dfchi_dv(i) = - cafor(i) * gca_high(i)
          dfchi_dchi(i) = - betchi(i)
         end do

       do i = 1, numcomp
c Note possible increase in rate at which AHP current develops
c       alpham_ahp(i) = dmin1(0.2d-4 * chi(i,L),0.01d0)
        alpham_ahp(i) = dmin1(1.0d-4 * chi(i,L),0.01d0)
        if (chi(i,L).le.500.d0) then
c         alpham_ahp_prime(i) = 0.2d-4
          alpham_ahp_prime(i) = 1.0d-4
        else
          alpham_ahp_prime(i) = 0.d0
        endif
       end do

       do i = 1, numcomp
        fmkahp(i) = alpham_ahp(i) * (1.d0 - mkahp(i,L))
c    x                  -.001d0 * mkahp(i,L)
     x                  -.010d0 * mkahp(i,L)
c       dfmkahp_dmkahp(i) = - alpham_ahp(i) - .001d0
        dfmkahp_dmkahp(i) = - alpham_ahp(i) - .010d0
        dfmkahp_dchi(i) = alpham_ahp_prime(i) *
     x                     (1.d0 - mkahp(i,L))
       end do

          do i = 1, numcomp

       K1 = IDNINT ( 4.d0 * (V(I,L) + 120.d0) )
       IF (K1.GT.640) K1 = 640
       IF (K1.LT.  0) K1 =   0

c      persistentNa_shift =  0.d0
c      persistentNa_shift =  8.d0
       persistentNa_shift = 10.d0
       K2 = IDNINT ( 4.d0 * (V(I,L)+persistentNa_shift+ 120.d0) )
       IF (K2.GT.640) K2 = 640
       IF (K2.LT.  0) K2 =   0

c            fastNa_shift = -2.0d0
c            fastNa_shift = -2.5d0
             fastNa_shift_SD = -3.5d0
             fastNa_shift_axon = fastNa_shift_SD + rel_axonshift 
       K0 = IDNINT ( 4.d0 * (V(I,L)+  fastNa_shift_SD+ 120.d0) )
       IF (K0.GT.640) K0 = 640
       IF (K0.LT.  0) K0 =   0
       K3 = IDNINT ( 4.d0 * (V(I,L)+  fastNa_shift_axon+ 120.d0) )
       IF (K3.GT.640) K3 = 640
       IF (K3.LT.  0) K3 =   0

         if (i.le.68) then   ! FOR SD
        fmnaf(i) = alpham_naf(k0) * (1.d0 - mnaf(i,L)) -
     X              betam_naf(k0) * mnaf(i,L)
        fhnaf(i) = alphah_naf(k0) * (1.d0 - hnaf(i,L)) -
     X              betah_naf(k0) * hnaf(i,L)
         else  ! for axon
        fmnaf(i) = alpham_naf(k3) * (1.d0 - mnaf(i,L)) -
     X              betam_naf(k3) * mnaf(i,L)
        fhnaf(i) = alphah_naf(k3) * (1.d0 - hnaf(i,L)) -
     X              betah_naf(k3) * hnaf(i,L)
         endif
        fmnap(i) = alpham_naf(k2) * (1.d0 - mnap(i,L)) -
     X              betam_naf(k2) * mnap(i,L)
        fmkdr(i) = alpham_kdr(k1) * (1.d0 - mkdr(i,L)) -
     X              betam_kdr(k1) * mkdr(i,L)
        fmka(i)  = alpham_ka (k1) * (1.d0 - mka(i,L)) -
     X              betam_ka (k1) * mka(i,L)
        fhka(i)  = alphah_ka (k1) * (1.d0 - hka(i,L)) -
     X              betah_ka (k1) * hka(i,L)
        fmk2(i)  = alpham_k2 (k1) * (1.d0 - mk2(i,L)) -
     X              betam_k2 (k1) * mk2(i,L)
        fhk2(i)  = alphah_k2 (k1) * (1.d0 - hk2(i,L)) -
     X              betah_k2 (k1) * hk2(i,L)
        fmkm(i)  = alpham_km (k1) * (1.d0 - mkm(i,L)) -
     X              betam_km (k1) * mkm(i,L)
        fmkc(i)  = alpham_kc (k1) * (1.d0 - mkc(i,L)) -
     X              betam_kc (k1) * mkc(i,L)
        fmcat(i) = alpham_cat(k1) * (1.d0 - mcat(i,L)) -
     X              betam_cat(k1) * mcat(i,L)
        fhcat(i) = alphah_cat(k1) * (1.d0 - hcat(i,L)) -
     X              betah_cat(k1) * hcat(i,L)
        fmcaL(i) = alpham_caL(k1) * (1.d0 - mcaL(i,L)) -
     X              betam_caL(k1) * mcaL(i,L)
        fmar(i)  = alpham_ar (k1) * (1.d0 - mar(i,L)) -
     X              betam_ar (k1) * mar(i,L)

       dfmnaf_dv(i) = dalpham_naf(k0) * (1.d0 - mnaf(i,L)) -
     X                  dbetam_naf(k0) * mnaf(i,L)
       dfmnap_dv(i) = dalpham_naf(k2) * (1.d0 - mnap(i,L)) -
     X                  dbetam_naf(k2) * mnap(i,L)
       dfhnaf_dv(i) = dalphah_naf(k1) * (1.d0 - hnaf(i,L)) -
     X                  dbetah_naf(k1) * hnaf(i,L)
       dfmkdr_dv(i) = dalpham_kdr(k1) * (1.d0 - mkdr(i,L)) -
     X                  dbetam_kdr(k1) * mkdr(i,L)
       dfmka_dv(i)  = dalpham_ka(k1) * (1.d0 - mka(i,L)) -
     X                  dbetam_ka(k1) * mka(i,L)
       dfhka_dv(i)  = dalphah_ka(k1) * (1.d0 - hka(i,L)) -
     X                  dbetah_ka(k1) * hka(i,L)
       dfmk2_dv(i)  = dalpham_k2(k1) * (1.d0 - mk2(i,L)) -
     X                  dbetam_k2(k1) * mk2(i,L)
       dfhk2_dv(i)  = dalphah_k2(k1) * (1.d0 - hk2(i,L)) -
     X                  dbetah_k2(k1) * hk2(i,L)
       dfmkm_dv(i)  = dalpham_km(k1) * (1.d0 - mkm(i,L)) -
     X                  dbetam_km(k1) * mkm(i,L)
       dfmkc_dv(i)  = dalpham_kc(k1) * (1.d0 - mkc(i,L)) -
     X                  dbetam_kc(k1) * mkc(i,L)
       dfmcat_dv(i) = dalpham_cat(k1) * (1.d0 - mcat(i,L)) -
     X                  dbetam_cat(k1) * mcat(i,L)
       dfhcat_dv(i) = dalphah_cat(k1) * (1.d0 - hcat(i,L)) -
     X                  dbetah_cat(k1) * hcat(i,L)
       dfmcaL_dv(i) = dalpham_caL(k1) * (1.d0 - mcaL(i,L)) -
     X                  dbetam_caL(k1) * mcaL(i,L)
       dfmar_dv(i)  = dalpham_ar(k1) * (1.d0 - mar(i,L)) -
     X                  dbetam_ar(k1) * mar(i,L)

       dfmnaf_dmnaf(i) =  - alpham_naf(k0) - betam_naf(k0)
       dfmnap_dmnap(i) =  - alpham_naf(k2) - betam_naf(k2)
       dfhnaf_dhnaf(i) =  - alphah_naf(k1) - betah_naf(k1)
       dfmkdr_dmkdr(i) =  - alpham_kdr(k1) - betam_kdr(k1)
       dfmka_dmka(i)  =   - alpham_ka (k1) - betam_ka (k1)
       dfhka_dhka(i)  =   - alphah_ka (k1) - betah_ka (k1)
       dfmk2_dmk2(i)  =   - alpham_k2 (k1) - betam_k2 (k1)
       dfhk2_dhk2(i)  =   - alphah_k2 (k1) - betah_k2 (k1)
       dfmkm_dmkm(i)  =   - alpham_km (k1) - betam_km (k1)
       dfmkc_dmkc(i)  =   - alpham_kc (k1) - betam_kc (k1)
       dfmcat_dmcat(i) =  - alpham_cat(k1) - betam_cat(k1)
       dfhcat_dhcat(i) =  - alphah_cat(k1) - betah_cat(k1)
       dfmcaL_dmcaL(i) =  - alpham_caL(k1) - betam_caL(k1)
       dfmar_dmar(i)  =   - alpham_ar (k1) - betam_ar (k1)

          end do

       dt2 = 0.5d0 * dt * dt

        do i = 1, numcomp
          v(i,L) = v(i,L) + dt * fv(i)
           do j = 1, numcomp
        v(i,L) = v(i,L) + dt2 * dfv_dv(i,j) * fv(j)
           end do
        v(i,L) = v(i,L) + dt2 * ( dfv_dchi(i) * fchi(i)
     X          + dfv_dmnaf(i) * fmnaf(i)
     X          + dfv_dmnap(i) * fmnap(i)
     X          + dfv_dhnaf(i) * fhnaf(i)
     X          + dfv_dmkdr(i) * fmkdr(i)
     X          + dfv_dmka(i)  * fmka(i)
     X          + dfv_dhka(i)  * fhka(i)
     X          + dfv_dmk2(i)  * fmk2(i)
     X          + dfv_dhk2(i)  * fhk2(i)
     X          + dfv_dmkm(i)  * fmkm(i)
     X          + dfv_dmkc(i)  * fmkc(i)
     X          + dfv_dmkahp(i)* fmkahp(i)
     X          + dfv_dmcat(i)  * fmcat(i)
     X          + dfv_dhcat(i) * fhcat(i)
     X          + dfv_dmcaL(i) * fmcaL(i)
     X          + dfv_dmar(i)  * fmar(i) )

        chi(i,L) = chi(i,L) + dt * fchi(i) + dt2 *
     X   (dfchi_dchi(i) * fchi(i) + dfchi_dv(i) * fv(i))
        mnaf(i,L) = mnaf(i,L) + dt * fmnaf(i) + dt2 *
     X   (dfmnaf_dmnaf(i) * fmnaf(i) + dfmnaf_dv(i)*fv(i))
        mnap(i,L) = mnap(i,L) + dt * fmnap(i) + dt2 *
     X   (dfmnap_dmnap(i) * fmnap(i) + dfmnap_dv(i)*fv(i))
        hnaf(i,L) = hnaf(i,L) + dt * fhnaf(i) + dt2 *
     X   (dfhnaf_dhnaf(i) * fhnaf(i) + dfhnaf_dv(i)*fv(i))
        mkdr(i,L) = mkdr(i,L) + dt * fmkdr(i) + dt2 *
     X   (dfmkdr_dmkdr(i) * fmkdr(i) + dfmkdr_dv(i)*fv(i))
        mka(i,L) =  mka(i,L) + dt * fmka(i) + dt2 *
     X   (dfmka_dmka(i) * fmka(i) + dfmka_dv(i) * fv(i))
        hka(i,L) =  hka(i,L) + dt * fhka(i) + dt2 *
     X   (dfhka_dhka(i) * fhka(i) + dfhka_dv(i) * fv(i))
        mk2(i,L) =  mk2(i,L) + dt * fmk2(i) + dt2 *
     X   (dfmk2_dmk2(i) * fmk2(i) + dfmk2_dv(i) * fv(i))
        hk2(i,L) =  hk2(i,L) + dt * fhk2(i) + dt2 *
     X   (dfhk2_dhk2(i) * fhk2(i) + dfhk2_dv(i) * fv(i))
        mkm(i,L) =  mkm(i,L) + dt * fmkm(i) + dt2 *
     X   (dfmkm_dmkm(i) * fmkm(i) + dfmkm_dv(i) * fv(i))
        mkc(i,L) =  mkc(i,L) + dt * fmkc(i) + dt2 *
     X   (dfmkc_dmkc(i) * fmkc(i) + dfmkc_dv(i) * fv(i))
        mkahp(i,L) = mkahp(i,L) + dt * fmkahp(i) + dt2 *
     X (dfmkahp_dmkahp(i)*fmkahp(i) + dfmkahp_dchi(i)*fchi(i))
        mcat(i,L) =  mcat(i,L) + dt * fmcat(i) + dt2 *
     X   (dfmcat_dmcat(i) * fmcat(i) + dfmcat_dv(i) * fv(i))
        hcat(i,L) =  hcat(i,L) + dt * fhcat(i) + dt2 *
     X   (dfhcat_dhcat(i) * fhcat(i) + dfhcat_dv(i) * fv(i))
        mcaL(i,L) =  mcaL(i,L) + dt * fmcaL(i) + dt2 *
     X   (dfmcaL_dmcaL(i) * fmcaL(i) + dfmcaL_dv(i) * fv(i))
        mar(i,L) =   mar(i,L) + dt * fmar(i) + dt2 *
     X   (dfmar_dmar(i) * fmar(i) + dfmar_dv(i) * fv(i))
c            endif
         end do

! Add membrane currents into membcurr for appropriate compartments
          do i = 1, 1 ! omit some basal comps
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do
c         do i = 14, 21
c          j = level(i)
c          membcurr(j) = membcurr(j) + fv(i) * c(i)
c         end do
c         do i = 26, 33
c          j = level(i)
c          membcurr(j) = membcurr(j) + fv(i) * c(i)
c         end do
          do i = 39, 68
           j = level(i)
           membcurr(j) = membcurr(j) + fv(i) * c(i)
          end do

            end do
c Finish loop L = 1 to numcell

         field_sup = 0.d0
         field_deep = 0.d0

         do i = 1, 12
        field_sup = field_sup + membcurr(i) / dabs(100.d0 - depth(i))
        field_deep = field_deep + membcurr(i) / dabs(500.d0 - depth(i))
         end do

2001          CONTINUE

6000    END



C  SETS UP TABLES FOR RATE FUNCTIONS
       SUBROUTINE SCORT_SETUP_semilunar
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)
      INTEGER I,J,K
      real*8 minf, hinf, taum, tauh, V, Z, shift_hnaf,
     X  shift_mkdr,
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
C FOR VOLTAGE, RANGE IS -120 TO +40 MV (absol.), 0.25 MV RESOLUTION


       DO 1, I = 0, 640
          V = dble(I)
          V = (V / 4.d0) - 120.d0

c gNa
           minf = 1.d0/(1.d0 + dexp((-V-38.d0)/10.d0))
           if (v.le.-30.d0) then
            taum = .025d0 + .14d0*dexp((v+30.d0)/10.d0)
           else
            taum = .02d0 + .145d0*dexp((-v-30.d0)/10.d0)
           endif
c from principal c. data, Martina & Jonas 1997, tau x 0.5
c Note that minf about the same for interneuron & princ. cell.
           alpham_naf(i) = minf / taum
           betam_naf(i) = 1.d0/taum - alpham_naf(i)

            shift_hnaf =  0.d0
        hinf = 1.d0/(1.d0 +
     x     dexp((v + shift_hnaf + 62.9d0)/10.7d0))
        tauh = 0.15d0 + 1.15d0/(1.d0+dexp((v+37.d0)/15.d0))
c from princ. cell data, Martina & Jonas 1997, tau x 0.5
            alphah_naf(i) = hinf / tauh
            betah_naf(i) = 1.d0/tauh - alphah_naf(i)

          shift_mkdr = 0.d0
c delayed rectifier, non-inactivating
       minf = 1.d0/(1.d0+dexp((-v-shift_mkdr-29.5d0)/10.0d0))
            if (v.le.-10.d0) then
             taum = .25d0 + 4.35d0*dexp((v+10.d0)/10.d0)
            else
             taum = .25d0 + 4.35d0*dexp((-v-10.d0)/10.d0)
            endif
              alpham_kdr(i) = minf / taum
              betam_kdr(i) = 1.d0 /taum - alpham_kdr(i)
c from Martina, Schultz et al., 1998. See espec. Table 1.

c A current: Huguenard & McCormick 1992, J Neurophysiol (TCR)
            minf = 1.d0/(1.d0 + dexp((-v-60.d0)/8.5d0))
            hinf = 1.d0/(1.d0 + dexp((v+78.d0)/6.d0))
        taum = .185d0 + .5d0/(dexp((v+35.8d0)/19.7d0) +
     x                            dexp((-v-79.7d0)/12.7d0))
        if (v.le.-63.d0) then
         tauh = .5d0/(dexp((v+46.d0)/5.d0) +
     x                  dexp((-v-238.d0)/37.5d0))
        else
         tauh = 9.5d0
        endif
           alpham_ka(i) = minf/taum
           betam_ka(i) = 1.d0 / taum - alpham_ka(i)
           alphah_ka(i) = hinf / tauh
           betah_ka(i) = 1.d0 / tauh - alphah_ka(i)

c h-current (anomalous rectifier), Huguenard & McCormick, 1992
           minf = 1.d0/(1.d0 + dexp((v+75.d0)/5.5d0))
           taum = 1.d0/(dexp(-14.6d0 -0.086d0*v) +
     x                   dexp(-1.87 + 0.07d0*v))
           alpham_ar(i) = minf / taum
           betam_ar(i) = 1.d0 / taum - alpham_ar(i)

c K2 K-current, McCormick & Huguenard
             minf = 1.d0/(1.d0 + dexp((-v-10.d0)/17.d0))
             hinf = 1.d0/(1.d0 + dexp((v+58.d0)/10.6d0))
            taum = 4.95d0 + 0.5d0/(dexp((v-81.d0)/25.6d0) +
     x                  dexp((-v-132.d0)/18.d0))
            tauh = 60.d0 + 0.5d0/(dexp((v-1.33d0)/200.d0) +
     x                  dexp((-v-130.d0)/7.1d0))
             alpham_k2(i) = minf / taum
             betam_k2(i) = 1.d0/taum - alpham_k2(i)
             alphah_k2(i) = hinf / tauh
             betah_k2(i) = 1.d0 / tauh - alphah_k2(i)

c voltage part of C-current, using 1994 kinetics, shift 60 mV
              if (v.le.-10.d0) then
       alpham_kc(i) = (2.d0/37.95d0)*dexp((v+50.d0)/11.d0 -
     x                                     (v+53.5)/27.d0)
       betam_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)-alpham_kc(i)
               else
       alpham_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)
       betam_kc(i) = 0.d0
               endif

c high-threshold gCa, from 1994, with 60 mV shift & no inactivn.
            alpham_cal(i) = 1.6d0/(1.d0+dexp(-.072d0*(v-5.d0)))
            betam_cal(i) = 0.1d0 * ((v+8.9d0)/5.d0) /
     x          (dexp((v+8.9d0)/5.d0) - 1.d0)

c M-current, from plast.f, with 60 mV shift
        alpham_km(i) = .02d0/(1.d0+dexp((-v-20.d0)/5.d0))
        betam_km(i) = .01d0 * dexp((-v-43.d0)/18.d0)

c T-current, from Destexhe, Neubig et al., 1998
         minf = 1.d0/(1.d0 + dexp((-v-56.d0)/6.2d0))
         hinf = 1.d0/(1.d0 + dexp((v+80.d0)/4.d0))
         taum = 0.204d0 + .333d0/(dexp((v+15.8d0)/18.2d0) +
     x                  dexp((-v-131.d0)/16.7d0))
          if (v.le.-81.d0) then
         tauh = 0.333 * dexp((v+466.d0)/66.6d0)
          else
         tauh = 9.32d0 + 0.333d0*dexp((-v-21.d0)/10.5d0)
          endif
              alpham_cat(i) = minf / taum
              betam_cat(i) = 1.d0/taum - alpham_cat(i)
              alphah_cat(i) = hinf / tauh
              betah_cat(i) = 1.d0 / tauh - alphah_cat(i)

1        CONTINUE

         do  i = 0, 639

      dalpham_naf(i) = (alpham_naf(i+1)-alpham_naf(i))/.25d0
      dbetam_naf(i) = (betam_naf(i+1)-betam_naf(i))/.25d0
      dalphah_naf(i) = (alphah_naf(i+1)-alphah_naf(i))/.25d0
      dbetah_naf(i) = (betah_naf(i+1)-betah_naf(i))/.25d0
      dalpham_kdr(i) = (alpham_kdr(i+1)-alpham_kdr(i))/.25d0
      dbetam_kdr(i) = (betam_kdr(i+1)-betam_kdr(i))/.25d0
      dalpham_ka(i) = (alpham_ka(i+1)-alpham_ka(i))/.25d0
      dbetam_ka(i) = (betam_ka(i+1)-betam_ka(i))/.25d0
      dalphah_ka(i) = (alphah_ka(i+1)-alphah_ka(i))/.25d0
      dbetah_ka(i) = (betah_ka(i+1)-betah_ka(i))/.25d0
      dalpham_k2(i) = (alpham_k2(i+1)-alpham_k2(i))/.25d0
      dbetam_k2(i) = (betam_k2(i+1)-betam_k2(i))/.25d0
      dalphah_k2(i) = (alphah_k2(i+1)-alphah_k2(i))/.25d0
      dbetah_k2(i) = (betah_k2(i+1)-betah_k2(i))/.25d0
      dalpham_km(i) = (alpham_km(i+1)-alpham_km(i))/.25d0
      dbetam_km(i) = (betam_km(i+1)-betam_km(i))/.25d0
      dalpham_kc(i) = (alpham_kc(i+1)-alpham_kc(i))/.25d0
      dbetam_kc(i) = (betam_kc(i+1)-betam_kc(i))/.25d0
      dalpham_cat(i) = (alpham_cat(i+1)-alpham_cat(i))/.25d0
      dbetam_cat(i) = (betam_cat(i+1)-betam_cat(i))/.25d0
      dalphah_cat(i) = (alphah_cat(i+1)-alphah_cat(i))/.25d0
      dbetah_cat(i) = (betah_cat(i+1)-betah_cat(i))/.25d0
      dalpham_caL(i) = (alpham_cal(i+1)-alpham_cal(i))/.25d0
      dbetam_caL(i) = (betam_cal(i+1)-betam_cal(i))/.25d0
      dalpham_ar(i) = (alpham_ar(i+1)-alpham_ar(i))/.25d0
      dbetam_ar(i) = (betam_ar(i+1)-betam_ar(i))/.25d0
       end do
2      CONTINUE

         do i = 640, 640
      dalpham_naf(i) =  dalpham_naf(i-1)
      dbetam_naf(i) =  dbetam_naf(i-1)
      dalphah_naf(i) = dalphah_naf(i-1)
      dbetah_naf(i) = dbetah_naf(i-1)
      dalpham_kdr(i) =  dalpham_kdr(i-1)
      dbetam_kdr(i) =  dbetam_kdr(i-1)
      dalpham_ka(i) =  dalpham_ka(i-1)
      dbetam_ka(i) =  dbetam_ka(i-1)
      dalphah_ka(i) =  dalphah_ka(i-1)
      dbetah_ka(i) =  dbetah_ka(i-1)
      dalpham_k2(i) =  dalpham_k2(i-1)
      dbetam_k2(i) =  dbetam_k2(i-1)
      dalphah_k2(i) =  dalphah_k2(i-1)
      dbetah_k2(i) =  dbetah_k2(i-1)
      dalpham_km(i) =  dalpham_km(i-1)
      dbetam_km(i) =  dbetam_km(i-1)
      dalpham_kc(i) =  dalpham_kc(i-1)
      dbetam_kc(i) =  dbetam_kc(i-1)
      dalpham_cat(i) =  dalpham_cat(i-1)
      dbetam_cat(i) =  dbetam_cat(i-1)
      dalphah_cat(i) =  dalphah_cat(i-1)
      dbetah_cat(i) =  dbetah_cat(i-1)
      dalpham_caL(i) =  dalpham_caL(i-1)
      dbetam_caL(i) =  dbetam_caL(i-1)
      dalpham_ar(i) =  dalpham_ar(i-1)
      dbetam_ar(i) =  dbetam_ar(i-1)
       end do   

4000   END

        SUBROUTINE SCORTMAJ_semilunar
C BRANCHED ACTIVE DENDRITES
     X             (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM,depth,level)
c Conductances: leak gL, coupling g, delayed rectifier gKDR, A gKA,
c C gKC, AHP gKAHP, K2 gK2, M gKM, low thresh Ca gCAT, high thresh
c gCAL, fast Na gNAF, persistent Na gNAP, h or anom. rectif. gAR.
c Note VAR = equil. potential for anomalous rectifier.
c Soma = comp. 1; 10 dendrites each with 13 compartments, 6-comp. axon
c Drop "glc"-like terms, just using "gl"-like
c CAFOR corresponds to "phi" in Traub et al., 1994
c Consistent set of units: nF, mV, ms, nA, microS

       INTEGER, PARAMETER:: numcomp = 74
! numcomp here must be compatible with numcomp_suppyrRS in calling prog.
        REAL*8 C(numcomp),GL(numcomp), GAM(0:numcomp, 0:numcomp)
        REAL*8 GNAF(numcomp),GCAT(numcomp), GKAHP(numcomp)
        REAL*8 GKDR(numcomp),GKA(numcomp),GKC(numcomp)
        REAL*8 GK2(numcomp),GNAP(numcomp),GAR(numcomp)
        REAL*8 GKM(numcomp), gcal(numcomp), CDENS
        REAL*8 JACOB(numcomp,numcomp),RI_SD,RI_AXON,RM_SD,RM_AXON
        INTEGER LEVEL(numcomp)
        REAL*8 GNAF_DENS(0:12), GCAT_DENS(0:12), GKDR_DENS(0:12)
        REAL*8 GKA_DENS(0:12), GKC_DENS(0:12), GKAHP_DENS(0:12)
        REAL*8 GCAL_DENS(0:12), GK2_DENS(0:12), GKM_DENS(0:12)
        REAL*8 GNAP_DENS(0:12), GAR_DENS(0:12)
        REAL*8 RES, RINPUT, Z, ELEN(numcomp)
        REAL*8 RSOMA, PI, BETCHI(numcomp), CAFOR(numcomp)
        REAL*8 RAD(numcomp), LEN(numcomp), GAM1, GAM2
        REAL*8 RIN, D(numcomp), AREA(numcomp), RI
        INTEGER NEIGH(numcomp,10), NNUM(numcomp), i, j, k, it
C FOR ESTABLISHING TOPOLOGY OF COMPARTMENTS
        real*8 depth(12) ! depth in microns of levels 1-12, assuming soma 
! at depth 500 microns 

        depth(1) = 300.d0
        depth(2) = 250.d0 ! now just obliques
        depth(3) = 250.d0 ! now just obliques
        depth(4) = 250.d0 ! now just obliques
        depth(5) = 250.d0
        depth(6) = 210.d0
        depth(7) = 170.d0
        depth(8) = 130.d0
        depth(9) =  90.d0
        depth(10) =  80.d0
        depth(11) =  70.d0
        depth(12) =  50.d0

        RI_SD = 250.d0
        RM_SD = 50000.d0
        RI_AXON = 100.d0
        RM_AXON = 1000.d0
        CDENS = 0.9d0

        PI = 3.14159d0

       do i = 0, 12
        gnaf_dens(i) = 10.d0
       end do
c       gnaf_dens(0) = 400.d0
!       gnaf_dens(0) = 120.d0
        gnaf_dens(0) = 200.d0
        gnaf_dens(1) = 120.d0
        gnaf_dens(2) =  75.d0
        gnaf_dens(5) = 100.d0
        gnaf_dens(6) =  75.d0

       do i = 0, 12
        gkdr_dens(i) = 0.d0
       end do
c       gkdr_dens(0) = 400.d0
c       gkdr_dens(0) = 100.d0
c       gkdr_dens(0) = 170.d0
        gkdr_dens(0) = 250.d0
c       gkdr_dens(1) = 100.d0
        gkdr_dens(1) = 150.d0
        gkdr_dens(2) =  75.d0
        gkdr_dens(5) = 100.d0
        gkdr_dens(6) =  75.d0

        gnap_dens(0) = 0.d0
        do i = 1, 12
          gnap_dens(i) = 0.0040d0 * gnaf_dens(i)
c         gnap_dens(i) = 0.002d0 * gnaf_dens(i)
c         gnap_dens(i) = 0.0030d0 * gnaf_dens(i)
        end do

        gcat_dens(0) = 0.d0
        do i = 1, 12
c         gcat_dens(i) = 0.5d0
          gcat_dens(i) = 0.1d0
        end do

        gcaL_dens(0) = 0.d0
        do i = 1, 6
          gcaL_dens(i) = 0.5d0
        end do
        do i = 7, 12
          gcaL_dens(i) = 0.5d0
        end do

       do i = 0, 12
        gka_dens(i) = 2.d0
       end do
        gka_dens(0) =100.d0 ! NOTE
        gka_dens(1) = 30.d0
        gka_dens(5) = 30.d0

      do i = 0, 12
c        gkc_dens(i)  = 12.00d0
         gkc_dens(i)  =  0.00d0
c        gkc_dens(i)  =  2.00d0
c        gkc_dens(i)  =  7.00d0
      end do
         gkc_dens(0) =  0.00d0
c        gkc_dens(1) = 7.5d0
c        gkc_dens(1) = 12.d0
         gkc_dens(1) = 15.d0
c        gkc_dens(2) = 7.5d0
         gkc_dens(2) = 10.d0
         gkc_dens(5) = 7.5d0
         gkc_dens(6) = 7.5d0

c       gkm_dens(0) = 2.d0 ! 9 Nov. 2005, see scort-pan.f of today
        gkm_dens(0) = 8.d0 ! 9 Nov. 2005, see scort-pan.f of today
! Above suppresses doublets, but still allows FRB with appropriate
! gNaP, gKC, and rel_axonshift (e.g. 6 mV)
        do i = 1, 12
         gkm_dens(i) = 2.5d0 * 1.50d0
        end do

        do i = 0, 12
c       gk2_dens(i) = 1.d0
        gk2_dens(i) = 0.1d0
        end do
        gk2_dens(0) = 0.d0

        gkahp_dens(0) = 0.d0
        do i = 1, 12
c        gkahp_dens(i) = 0.200d0
         gkahp_dens(i) = 0.100d0
c        gkahp_dens(i) = 0.050d0
        end do

        gar_dens(0) = 0.d0
        do i = 1, 12
         gar_dens(i) = 0.25d0
        end do

c       WRITE   (6,9988)
9988    FORMAT(2X,'I',4X,'NADENS',' CADENS(T)',' KDRDEN',' KAHPDE',
     X     ' KCDENS',' KADENS')
        DO 9989, I = 0, 12
c         WRITE (6,9990) I, gnaf_dens(i), gcat_dens(i), gkdr_dens(i),
c    X  gkahp_dens(i), gkc_dens(i), gka_dens(i)
9990    FORMAT(2X,I2,2X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2)
9989    CONTINUE


        level(1) = 1
        do i = 2, 13
         level(i) = 2
        end do
        do i = 14, 25
           level(i) = 3
        end do
        do i = 26, 37
           level(i) = 4
        end do
        level(38) = 5
        level(39) = 6
        level(40) = 7
        level(41) = 8
        level(42) = 8
        level(43) = 9
        level(44) = 9
        do i = 45, 52
           level(i) = 10
        end do
        do i = 53, 60
           level(i) = 11
        end do
        do i = 61, 68
           level(i) = 12
        end do

        do i =  69, 74
         level(i) = 0
        end do

c connectivity of axon
        nnum( 69) = 2
        nnum( 70) = 3
        nnum( 71) = 3
        nnum( 73) = 3
        nnum( 72) = 1
        nnum( 74) = 1
         neigh(69,1) =  1
         neigh(69,2) = 70
         neigh(70,1) = 69
         neigh(70,2) = 71
         neigh(70,3) = 73
         neigh(71,1) = 70
         neigh(71,2) = 72
         neigh(71,3) = 73
         neigh(73,1) = 70
         neigh(73,2) = 71
         neigh(73,3) = 74
         neigh(72,1) = 71
         neigh(74,1) = 73

c connectivity of SD part
          nnum(1) = 10
          neigh(1,1) = 69
          neigh(1,2) =  2
          neigh(1,3) =  3
          neigh(1,4) =  4
          neigh(1,5) =  5
          neigh(1,6) =  6
          neigh(1,7) =  7
          neigh(1,8) =  8
          neigh(1,9) =  9
          neigh(1,10) = 38

          do i = 2, 9
           nnum(i) = 2
           neigh(i,1) = 1
           neigh(i,2) = i + 12
          end do

          do i = 14, 21
            nnum(i) = 2
            neigh(i,1) = i - 12
            neigh(i,2) = i + 12
          end do

          do i = 26, 33
            nnum(i) = 1
            neigh(i,1) = i - 12
          end do

          do i = 10, 13
            nnum(i) = 2
            neigh(i,1) = 38
            neigh(i,2) = i + 12
          end do

          do i = 22, 25
            nnum(i) = 2
            neigh(i,1) = i - 12
            neigh(i,2) = i + 12
          end do

          do i = 34, 37
            nnum(i) = 1
            neigh(i,1) = i - 12
          end do

          nnum(38) = 6
          neigh(38,1) = 1
          neigh(38,2) = 39
          neigh(38,3) = 10
          neigh(38,4) = 11
          neigh(38,5) = 12
          neigh(38,6) = 13

          nnum(39) = 2
          neigh(39,1) = 38
          neigh(39,2) = 40

          nnum(40) = 3
          neigh(40,1) = 39
          neigh(40,2) = 41
          neigh(40,3) = 42

          nnum(41) = 3
          neigh(41,1) = 40
          neigh(41,2) = 42
          neigh(41,3) = 43

          nnum(42) = 3
          neigh(42,1) = 40
          neigh(42,2) = 41
          neigh(42,3) = 44

           nnum(43) = 5
           neigh(43,1) = 41
           neigh(43,2) = 45
           neigh(43,3) = 46
           neigh(43,4) = 47
           neigh(43,5) = 48

           nnum(44) = 5
           neigh(44,1) = 42
           neigh(44,2) = 49
           neigh(44,3) = 50
           neigh(44,4) = 51
           neigh(44,5) = 52

           nnum(45) = 5
           neigh(45,1) = 43
           neigh(45,2) = 53
           neigh(45,3) = 46
           neigh(45,4) = 47
           neigh(45,5) = 48

           nnum(46) = 5
           neigh(46,1) = 43
           neigh(46,2) = 54
           neigh(46,3) = 45
           neigh(46,4) = 47
           neigh(46,5) = 48

           nnum(47) = 5
           neigh(47,1) = 43
           neigh(47,2) = 55
           neigh(47,3) = 45
           neigh(47,4) = 46
           neigh(47,5) = 48

           nnum(48) = 5
           neigh(48,1) = 43
           neigh(48,2) = 56
           neigh(48,3) = 45
           neigh(48,4) = 46
           neigh(48,5) = 47

           nnum(49) = 5
           neigh(49,1) = 44
           neigh(49,2) = 57
           neigh(49,3) = 50
           neigh(49,4) = 51
           neigh(49,5) = 52

           nnum(50) = 5
           neigh(50,1) = 44
           neigh(50,2) = 58
           neigh(50,3) = 49
           neigh(50,4) = 51
           neigh(50,5) = 52

           nnum(51) = 5
           neigh(51,1) = 44
           neigh(51,2) = 59
           neigh(51,3) = 49
           neigh(51,4) = 50
           neigh(51,5) = 52

           nnum(52) = 5
           neigh(52,1) = 44
           neigh(52,2) = 60
           neigh(52,3) = 49
           neigh(52,4) = 51
           neigh(52,5) = 50

          do i = 53, 60
           nnum(i) = 2
           neigh(i,1) = i - 8
           neigh(i,2) = i + 8
          end do

          do i = 61, 68
           nnum(i) = 1
           neigh(i,1) = i - 8
          end do

c        DO 332, I = 1, 74
         DO I = 1, 74
c          WRITE(6,3330) I, NEIGH(I,1),NEIGH(I,2),NEIGH(I,3),NEIGH(I,4),
c    X NEIGH(I,5),NEIGH(I,6),NEIGH(I,7),NEIGH(I,8),NEIGH(I,9),
c    X NEIGH(I,10)
3330     FORMAT(2X,11I5)
         END DO
332      CONTINUE
c         DO 858, I = 1, 74
          DO I = 1, 74
c          DO 858, J = 1, NNUM(I)
           DO J = 1, NNUM(I)
            K = NEIGH(I,J)
            IT = 0
c           DO 859, L = 1, NNUM(K)
            DO  L = 1, NNUM(K)
             IF (NEIGH(K,L).EQ.I) IT = 1
            END DO
859         CONTINUE
             IF (IT.EQ.0) THEN
c             WRITE(6,8591) I, K
8591          FORMAT(' ASYMMETRY IN NEIGH MATRIX ',I4,I4)
              STOP
             ENDIF
          END DO
          END DO
858       CONTINUE

c length and radius of axonal compartments
c Note shortened "initial segment"
          len(69) = 25.d0
          do i = 70, 74
            len(i) = 50.d0
          end do
          rad( 69) = 0.90d0
c         rad( 69) = 0.80d0
          rad( 70) = 0.7d0
          do i = 71, 74
           rad(i) = 0.5d0
          end do

c  length and radius of SD compartments
          len(1) = 15.d0
          rad(1) =  8.d0

          do i = 2, 68
           len(i) = 50.d0
          end do

          do i = 2, 37
            rad(i) = 0.5d0
          end do

          z = 4.0d0
          rad(38) = z
          rad(39) = 0.9d0 * z
          rad(40) = 0.8d0 * z
          rad(41) = 0.5d0 * z
          rad(42) = 0.5d0 * z
          rad(43) = 0.5d0 * z
          rad(44) = 0.5d0 * z
          do i = 45, 68
           rad(i) = 0.2d0 * z
          end do


c       WRITE(6,919)
919     FORMAT('COMPART.',' LEVEL ',' RADIUS ',' LENGTH(MU)')
c       DO 920, I = 1, 74
c920      WRITE(6,921) I, LEVEL(I), RAD(I), LEN(I)
921     FORMAT(I3,5X,I2,3X,F6.2,1X,F6.1,2X,F4.3)

        DO 120, I = 1, 74
          AREA(I) = 2.d0 * PI * RAD(I) * LEN(I)
      if((i.gt.1).and.(i.le.68)) area(i) = 2.d0 * area(i)
C    CORRECTION FOR CONTRIBUTION OF SPINES TO AREA
          K = LEVEL(I)
          C(I) = CDENS * AREA(I) * (1.D-8)

           if (k.ge.1) then
          GL(I) = (1.D-2) * AREA(I) / RM_SD
           else
          GL(I) = (1.D-2) * AREA(I) / RM_AXON
           endif

          GNAF(I) = GNAF_DENS(K) * AREA(I) * (1.D-5)
          GNAP(I) = GNAP_DENS(K) * AREA(I) * (1.D-5)
          GCAT(I) = GCAT_DENS(K) * AREA(I) * (1.D-5)
          GKDR(I) = GKDR_DENS(K) * AREA(I) * (1.D-5)
          GKA(I) = GKA_DENS(K) * AREA(I) * (1.D-5)
          GKC(I) = GKC_DENS(K) * AREA(I) * (1.D-5)
          GKAHP(I) = GKAHP_DENS(K) * AREA(I) * (1.D-5)
          GKM(I) = GKM_DENS(K) * AREA(I) * (1.D-5)
          GCAL(I) = GCAL_DENS(K) * AREA(I) * (1.D-5)
          GK2(I) = GK2_DENS(K) * AREA(I) * (1.D-5)
          GAR(I) = GAR_DENS(K) * AREA(I) * (1.D-5)
c above conductances should be in microS
120           continue

         Z = 0.d0
c        DO 1019, I = 2, 68
         DO I = 2, 68
           Z = Z + AREA(I)
         END DO
1019     CONTINUE
c        WRITE(6,1020) Z
1020     FORMAT(2X,' TOTAL DENDRITIC AREA ',F7.0)

c       DO 140, I = 1, 74
        DO I = 1, 74
c       DO 140, K = 1, NNUM(I)
        DO K = 1, NNUM(I)
         J = NEIGH(I,K)
           if (level(i).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM1 =100.d0 * PI * RAD(I) * RAD(I) / ( RI * LEN(I) )

           if (level(j).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM2 =100.d0 * PI * RAD(J) * RAD(J) / ( RI * LEN(J) )
         GAM(I,J) = 2.d0/( (1.d0/GAM1) + (1.d0/GAM2) )
	 END DO
	 END DO
c DISCONNECT BASAL DENDRITES FROM SOMA
         do i = 2, 9
          gam(1,i) = 0.d0
          gam(i,1) = 0.d0
         end do

140     CONTINUE
c gam computed in microS

c       DO 299, I = 1, 74
        DO I = 1, 74
299       BETCHI(I) = .05d0
        END DO
        BETCHI( 1) =  .01d0

c       DO 300, I = 1, 74
        DO I = 1, 74
c300     D(I) = 2.D-4
300     D(I) = 5.D-4
        END DO
c       DO 301, I = 1, 74
        DO I = 1, 74
         IF (LEVEL(I).EQ.1) D(I) = 2.D-3
        END DO
301     CONTINUE
C  NOTE NOTE NOTE  (DIFFERENT FROM SWONG)


c      DO 160, I = 1, 74
       DO I = 1, 74
160     CAFOR(I) = 5200.d0 / (AREA(I) * D(I))
       END DO
C     NOTE CORRECTION

c       do 200, i = 1, 74
        do i = 1, numcomp
200     C(I) = 1000.d0 * C(I)
        end do
C     TO GO FROM MICROF TO NF.

c     DO 909, I = 1, 74
      DO I = 1, numcomp
       JACOB(I,I) = - GL(I)
c     DO 909, J = 1, NNUM(I)
      DO J = 1, NNUM(I)
         K = NEIGH(I,J)
         IF (I.EQ.K) THEN
c            WRITE(6,510) I
510          FORMAT(' UNEXPECTED SYMMETRY IN NEIGH ',I4)
         ENDIF
         JACOB(I,K) = GAM(I,K)
         JACOB(I,I) = JACOB(I,I) - GAM(I,K)
       END DO
       END DO
909   CONTINUE

c 15 Jan. 2001: make correction for c(i)
          do i = 1, numcomp
          do j = 1, numcomp
             jacob(i,j) = jacob(i,j) / c(i)
          end do
          end do

c      DO 500, I = 1, 74
       DO I = 1, 74
c       WRITE (6,501) I,C(I)
501     FORMAT(1X,I3,' C(I) = ',F7.4)
       END DO
500     CONTINUE
        END


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

c      value = value * (10.2d0 * dexp(-(t-10.d0)/122.d0) +
c    &    1.1d0 * dexp(-(t-10.d0)/587.d0))

c      value = value * (10.2d0 * dexp(-(t-10.d0)/250.d0) +
       value = value * (10.2d0 * dexp(-(t-10.d0)/200.d0) +
     &    0.0d0 * dexp(-(t-10.d0)/587.d0))
              endif

                 end


! 21 March 2021, from piriform integrate_deepbask.f
! substitute scort_setup_suppyrRS.f and modify various
! parameters per results in /multipolar
c 23 Aug 2019, taken from son_of_groucho, for use in piriform.f
c Was deepbaskx.f, now renamed deepbask.f
! Integration program for superior & deep basket & axo-axonic cells
! From baskn.f in supergj.f

       SUBROUTINE integrate_multipolar (O, time, numcell, V, curr,
     &  initialize, firstcell, lastcell,
     &  gAMPA, gNMDA, gGABA_A, Mg, gapcon, totaxgj, gjtable, dt,
     &  chi,mnaf,mnap,
     &  hnaf,mkdr,mka,
     &  hka,mk2,hk2,
     &  mkm,mkc,mkahp,
     &  mcat,hcat,mcal,
     &  mar)

           SAVE

       integer, parameter:: numcomp = 59  ! should be compat. with calling prog

       integer numcell, totaxgj, gjtable(totaxgj,4)
       integer initialize, firstcell, lastcell
       INTEGER J1, I, J, K, L, L1, O, K1
       REAL*8  Z, Z1, Z2, curr(numcomp,numcell), c(numcomp)
       REAL*8 dt, time, Mg, gapcon
c Usual dt in this program .002 ms

c CINV is 1/C, i.e. inverse capacitance

       real*8 v(numcomp,numcell), chi(numcomp,numcell), cinv(numcomp),
     x mnaf(numcomp,numcell),hnaf(numcomp,numcell), 
     x mkdr(numcomp,numcell),
     x mka(numcomp,numcell),hka(numcomp,numcell),mk2(numcomp,numcell),
     x hk2(numcomp,numcell),mkm(numcomp,numcell),
     x mkc(numcomp,numcell),mkahp(numcomp,numcell),
     x mcat(numcomp,numcell),hcat(numcomp,numcell),
     x mcal(numcomp,numcell),mar(numcomp,numcell),
     x jacob(numcomp,numcomp),betchi(numcomp),
     x gam(0:numcomp,0:numcomp),gL(numcomp),gnaf(numcomp),
     x gnap(numcomp),gkdr(numcomp),gka(numcomp),
     x gk2(numcomp),gkm(numcomp),gkc(numcomp),gkahp(numcomp),
     x gcat(numcomp),gcaL(numcomp),gar(numcomp),
     x cafor(numcomp), ggaba_a(numcomp,numcell),
     x gampa(numcomp,numcell),gnmda(numcomp,numcell)
       real*8
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
       real*8 vL,vk,vna,var,vca,vgaba_a

        INTEGER NEIGH(numcomp,5), NNUM(numcomp)
        real*8 fastna_shift

c the f's are the functions giving 1st derivatives for evolution of
c the differential equations for the voltages (v), calcium (chi), and
c other state variables.
       real*8 fv(numcomp), fchi(numcomp),fmnaf(numcomp),
     x fhnaf(numcomp),fmkdr(numcomp),
     x fmka(numcomp),fhka(numcomp),fmk2(numcomp),fhk2(numcomp),
     x fmkm(numcomp),fmkc(numcomp),fmkahp(numcomp),
     x fmcat(numcomp),fhcat(numcomp),fmcal(numcomp),fmar(numcomp)

c below are for calculating the partial derivatives
       real*8 dfv_dv(numcomp,numcomp), dfv_dchi(numcomp), 
     x  dfv_dmnaf(numcomp),
     x  dfv_dhnaf(numcomp),dfv_dmkdr(numcomp),
     x  dfv_dmka(numcomp),dfv_dhka(numcomp),
     x  dfv_dmk2(numcomp),dfv_dhk2(numcomp),
     x  dfv_dmkm(numcomp),dfv_dmkc(numcomp),
     x  dfv_dmkahp(numcomp),dfv_dmcat(numcomp),
     x  dfv_dhcat(numcomp),dfv_dmcal(numcomp),
     x  dfv_dmar(numcomp)

        real*8 dfchi_dv(numcomp), dfchi_dchi(numcomp),
     x dfmnaf_dmnaf(numcomp), dfmnaf_dv(numcomp),dfhnaf_dhnaf(numcomp),
     x dfhnaf_dv(numcomp),dfmkdr_dmkdr(numcomp),dfmkdr_dv(numcomp),
     x dfmka_dmka(numcomp),dfmka_dv(numcomp),
     x dfhka_dhka(numcomp),dfhka_dv(numcomp),
     x dfmk2_dmk2(numcomp),dfmk2_dv(numcomp),
     x dfhk2_dhk2(numcomp),dfhk2_dv(numcomp),
     x dfmkm_dmkm(numcomp),dfmkm_dv(numcomp),
     x dfmkc_dmkc(numcomp),dfmkc_dv(numcomp),
     x dfmcat_dmcat(numcomp),dfmcat_dv(numcomp),dfhcat_dhcat(numcomp),
     x dfhcat_dv(numcomp),dfmcal_dmcal(numcomp),dfmcal_dv(numcomp),
     x dfmar_dmar(numcomp),dfmar_dv(numcomp),dfmkahp_dchi(numcomp),
     x dfmkahp_dmkahp(numcomp), dt2

         INTEGER  K0
       REAL*8 OPEN(numcomp),gamma(numcomp),gamma_prime(numcomp)
c gamma is function of chi used in calculating KC conductance
       REAL*8 alpham_ahp(numcomp), alpham_ahp_prime(numcomp)
       REAL*8 gna_tot(numcomp),gk_tot(numcomp),gca_tot(numcomp)
       REAL*8 gca_high(numcomp), gar_tot(numcomp)
c this will be gCa conductance corresponding to high-thresh channels
       REAL*8 A, BB1, BB2

c do initialization on 1st time step
c      if (O.eq.1) then
       if (initialize.eq.0) then

c Program fnmda assumes A, BB1, BB2 defined in calling program
c as follows:
         A = DEXP(-2.847d0)
         BB1 = DEXP(-.693d0)
         BB2 = DEXP(-3.101d0)

c      CALL  DEEPBASK_SETUP
       CALL scort_setup_suppyrRS
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)

        CALL multipolarMAJ (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM)

          do i = 1, 59
             cinv(i) = 1.d0 / c(i)
          end do

C  IN MILLIMOLAR

        VL = -65.d0
        VK =  -85.d0
        VNA = 50.d0
        VCA = 125.d0
        VAR = -40.d0
        VGABA_A = -75.d0


c ? initialize membrane state variables?
        do L = 1, numcell
        do i = 1, numcomp
          v(i,L) = VL
	  chi(i,L) = 0.d0
	mnaf(i,L) = 0.d0
	mkdr(i,L) = 0.d0
	mk2(i,L) = 0.d0
	mkm(i,L) = 0.d0
	mkc(i,L) = 0.d0
	mkahp(i,L) = 0.d0
	mcat(i,L) = 0.d0
	mcal(i,L) = 0.d0
	mar(i,L) = 0.d0

        k1 = idnint (4.d0 * (vL + 120.d0))

      hnaf(i,L) = alphah_naf(k1)/(alphah_naf(k1)+betah_naf(k1))
      hka(i,L) = alphah_ka(k1)/(alphah_ka(k1)+betah_ka(k1))
      hk2(i,L) = alphah_k2(k1)/(alphah_k2(k1)+betah_k2(k1))
      hcat(i,L)=alphah_cat(k1)/(alphah_cat(k1)+betah_cat(k1))
         end do
         end do


          do i = 1, numcomp
c                 gnap(i) = 0.d0
                  gk2(i) = 0.d0
c                 gkm(i) = 0.d0
c                 gkahp(i) = 0.d0
c                 gcat(i) = 0.d0
c                 gar(i) = 0.d0
		  open(i) = 0.d0
          end do

              goto 1000

c End initialization
             endif

c          do L = 1, numcell
           do L = firstcell, lastcell


       DO I = 1, numcomp
          FV(I) = -GL(I) * (V(I,L) - VL) * cinv(i)
c         DO 302, J = 1, NNUM(I)
          DO J = 1, NNUM(I)
             K = NEIGH(I,J)
302     FV(I) = FV(I) + GAM(I,K) * (V(K,L) - V(I,L)) * cinv(i)
          END DO
        END DO
301    CONTINUE


        CALL FNMDA (V, OPEN, numcell, numcomp, MG, L, 
     &    A, BB1, BB2)

      DO i = 1, numcomp
421    FV(I) = FV(I) + ( CURR(I,L)
     X   - (gampa(I,L) + open(i) * gnmda(I,L))*V(I,L)
     X   - ggaba_a(I,L)*(V(I,L)-Vgaba_a) ) * cinv(i)
      END DO
c above assumes equil. potential for AMPA & NMDA = 0 mV

       do m = 1, totaxgj
        if (gjtable(m,1).eq.L) then
         L1 = gjtable(m,3)
         igap1 = gjtable(m,2)
         igap2 = gjtable(m,4)
 	fv(igap1) = fv(igap1) + gapcon *
     &   (v(igap2,L1) - v(igap1,L)) * cinv(igap1)
        else if (gjtable(m,3).eq.L) then
         L1 = gjtable(m,1)
         igap1 = gjtable(m,4)
         igap2 = gjtable(m,2)
 	fv(igap1) = fv(igap1) + gapcon *
     &   (v(igap2,L1) - v(igap1,L)) * cinv(igap1)
        endif
       end do ! do m

c      do i = 1, ngap_FS(L) ! obsolete gj code
c      L1 = list_gap_FS(L,i)
c       fv(dendsite) = fv(dendsite) + gapconid_FS *
c    &   (vdgap_global_FS(L1) - v(dendsite,L)) * cinv(dendsite)
c      end do  ! obsolete gj code

       do i = 1, numcomp
        gamma(i) = dmin1 (1.d0, .004d0 * chi(i,L))
        if (chi(i,L).le.250.d0) then
          gamma_prime(i) = .004d0
        else
          gamma_prime(i) = 0.d0
        endif
       end do

c     DO 88, I = 1, numcomp
      DO I = 1, numcomp
       gna_tot(i) = gnaf(i) * (mnaf(i,L)**3) * hnaf(i,L) +
     x     gnap(i) * (mnaf(i,L)**3)
       gk_tot(i) = gkdr(i) * (mkdr(i,L)**4) +
     x             gka(i)  * (mka(i,L)**4) * hka(i,L) +
     x             gk2(i)  * mk2(i,L) * hk2(i,L) +
     x             gkm(i)  * mkm(i,L) +
     x             gkc(i)  * mkc(i,L) * gamma(i) +
     x             gkahp(i)* mkahp(i,L)
       gca_tot(i) = gcat(i) * (mcat(i,L)**2) * hcat(i,L) +
     x              gcaL(i) * (mcaL(i,L)**2)
       gca_high(i) =
     x              gcaL(i) * (mcaL(i,L)**2)
       gar_tot(i) = gar(i) * mar(i,L)


88     FV(I) = FV(I) - ( gna_tot(i) * (v(i,L) - vna)
     X  + gk_tot(i) * (v(i,L) - vK)
     X  + gca_tot(i) * (v(i,L) - vCa)
     X  + gar_tot(i) * (v(i,L) - var) ) * cinv(i)
       END DO

         do i = 1, numcomp
         do j = 1, numcomp
          if (i.ne.j) then
            dfv_dv(i,j) = jacob(i,j)
          else
            dfv_dv(i,j) = jacob(i,i) - cinv(i) *
     X  (gna_tot(i) + gk_tot(i) + gca_tot(i) + gar_tot(i)
     X   + ggaba_a(i,L) + gampa(i,L)
     X   + open(i) * gnmda(I,L) )
          endif
         end do
         end do

          do i = 1, numcomp
        dfv_dchi(i)  = - cinv(i) * gkc(i) * mkc(i,L) *
     x                     gamma_prime(i) * (v(i,L)-vK)
        dfv_dmnaf(i) = -3.d0 * cinv(i) * (mnaf(i,L)**2) *
     X    (gnaf(i) * hnaf(i,L) + gnap(i)) * (v(i,L) - vna)
        dfv_dhnaf(i) = - cinv(i) * gnaf(i) * (mnaf(i,L)**3) *
     X                    (v(i,L) - vna)
        dfv_dmkdr(i) = -4.d0 * cinv(i)*gkdr(i) * (mkdr(i,L)**3)
     X                   * (v(i,L) - vK)
        dfv_dmka(i)  = -4.d0 * cinv(i)*gka(i) * (mka(i,L)**3) *
     X                   hka(i,L) * (v(i,L) - vK)
        dfv_dhka(i)  = - cinv(i) * gka(i) * (mka(i,L)**4) *
     X                    (v(i,L) - vK)
       dfv_dmk2(i)  = - cinv(i)*gk2(i) * hk2(i,L) * (v(i,L)-vK)
       dfv_dhk2(i)  = - cinv(i)*gk2(i) * mk2(i,L) * (v(i,L)-vK)
       dfv_dmkm(i)  = - cinv(i)*gkm(i) * (v(i,L) - vK)
       dfv_dmkc(i)  = - cinv(i)*gkc(i) * gamma(i) * (v(i,L)-vK)
       dfv_dmkahp(i)= - cinv(i)*gkahp(i) * (v(i,L) - vK)
       dfv_dmcat(i)  = -2.d0 * cinv(i) * gcat(i) * mcat(i,L) *
     X                    hcat(i,L) * (v(i,L) - vCa)
        dfv_dhcat(i) = - cinv(i) * gcat(i) * (mcat(i,L)**2) *
     X                  (v(i,L) - vCa)
        dfv_dmcal(i) = -2.d0 * cinv(i) * gcal(i) * mcal(i,L) *
     X                      (v(i,L) - vCa)
        dfv_dmar(i) = - cinv(i) * gar(i) * (v(i,L) - var)
          end do

         do i = 1, numcomp
          fchi(i) = - cafor(i) * gca_high(i) * (v(i,L) - vca)
     x       - betchi(i) * chi(i,L)
          dfchi_dv(i) = - cafor(i) * gca_high(i)
          dfchi_dchi(i) = - betchi(i)
         end do

       do i = 1, numcomp
        alpham_ahp(i) = dmin1(0.2d-4 * chi(i,L),0.01d0)
        if (chi(i,L).le.500.d0) then
          alpham_ahp_prime(i) = 0.2d-4
        else
          alpham_ahp_prime(i) = 0.d0
        endif
       end do

       do i = 1, numcomp
        fmkahp(i) = alpham_ahp(i) * (1.d0 - mkahp(i,L))
     x                  -.001d0 * mkahp(i,L)
        dfmkahp_dmkahp(i) = - alpham_ahp(i) - .001d0
        dfmkahp_dchi(i) = alpham_ahp_prime(i) *
     x                     (1.d0 - mkahp(i,L))
       end do

          do i = 1, numcomp


       K1 = IDNINT ( 4.d0 * (V(I,L) + 120.d0) )
       IF (K1.GT.640) K1 = 640
       IF (K1.LT.  0) K1 =   0

             fastNa_shift = -2.5d0
       K0 = IDNINT ( 4.d0 * (V(I,L)+  fastNa_shift+ 120.d0) )
       IF (K0.GT.640) K0 = 640
       IF (K0.LT.  0) K0 =   0


        fmnaf(i) = alpham_naf(k0) * (1.d0 - mnaf(i,L)) -
     X              betam_naf(k0) * mnaf(i,L)
        fhnaf(i) = alphah_naf(k1) * (1.d0 - hnaf(i,L)) -
     X              betah_naf(k1) * hnaf(i,L)
        fmkdr(i) = alpham_kdr(k1) * (1.d0 - mkdr(i,L)) -
     X              betam_kdr(k1) * mkdr(i,L)
        fmka(i)  = alpham_ka (k1) * (1.d0 - mka(i,L)) -
     X              betam_ka (k1) * mka(i,L)
        fhka(i)  = alphah_ka (k1) * (1.d0 - hka(i,L)) -
     X              betah_ka (k1) * hka(i,L)
        fmk2(i)  = alpham_k2 (k1) * (1.d0 - mk2(i,L)) -
     X              betam_k2 (k1) * mk2(i,L)
        fhk2(i)  = alphah_k2 (k1) * (1.d0 - hk2(i,L)) -
     X              betah_k2 (k1) * hk2(i,L)
        fmkm(i)  = alpham_km (k1) * (1.d0 - mkm(i,L)) -
     X              betam_km (k1) * mkm(i,L)
        fmkc(i)  = alpham_kc (k1) * (1.d0 - mkc(i,L)) -
     X              betam_kc (k1) * mkc(i,L)
        fmcat(i) = alpham_cat(k1) * (1.d0 - mcat(i,L)) -
     X              betam_cat(k1) * mcat(i,L)
        fhcat(i) = alphah_cat(k1) * (1.d0 - hcat(i,L)) -
     X              betah_cat(k1) * hcat(i,L)
        fmcaL(i) = alpham_caL(k1) * (1.d0 - mcaL(i,L)) -
     X              betam_caL(k1) * mcaL(i,L)
        fmar(i)  = alpham_ar (k1) * (1.d0 - mar(i,L)) -
     X              betam_ar (k1) * mar(i,L)

       dfmnaf_dv(i) = dalpham_naf(k0) * (1.d0 - mnaf(i,L)) -
     X                  dbetam_naf(k0) * mnaf(i,L)
       dfhnaf_dv(i) = dalphah_naf(k1) * (1.d0 - hnaf(i,L)) -
     X                  dbetah_naf(k1) * hnaf(i,L)
       dfmkdr_dv(i) = dalpham_kdr(k1) * (1.d0 - mkdr(i,L)) -
     X                  dbetam_kdr(k1) * mkdr(i,L)
       dfmka_dv(i)  = dalpham_ka(k1) * (1.d0 - mka(i,L)) -
     X                  dbetam_ka(k1) * mka(i,L)
       dfhka_dv(i)  = dalphah_ka(k1) * (1.d0 - hka(i,L)) -
     X                  dbetah_ka(k1) * hka(i,L)
       dfmk2_dv(i)  = dalpham_k2(k1) * (1.d0 - mk2(i,L)) -
     X                  dbetam_k2(k1) * mk2(i,L)
       dfhk2_dv(i)  = dalphah_k2(k1) * (1.d0 - hk2(i,L)) -
     X                  dbetah_k2(k1) * hk2(i,L)
       dfmkm_dv(i)  = dalpham_km(k1) * (1.d0 - mkm(i,L)) -
     X                  dbetam_km(k1) * mkm(i,L)
       dfmkc_dv(i)  = dalpham_kc(k1) * (1.d0 - mkc(i,L)) -
     X                  dbetam_kc(k1) * mkc(i,L)
       dfmcat_dv(i) = dalpham_cat(k1) * (1.d0 - mcat(i,L)) -
     X                  dbetam_cat(k1) * mcat(i,L)
       dfhcat_dv(i) = dalphah_cat(k1) * (1.d0 - hcat(i,L)) -
     X                  dbetah_cat(k1) * hcat(i,L)
       dfmcaL_dv(i) = dalpham_caL(k1) * (1.d0 - mcaL(i,L)) -
     X                  dbetam_caL(k1) * mcaL(i,L)
       dfmar_dv(i)  = dalpham_ar(k1) * (1.d0 - mar(i,L)) -
     X                  dbetam_ar(k1) * mar(i,L)

       dfmnaf_dmnaf(i) =  - alpham_naf(k0) - betam_naf(k0)
       dfhnaf_dhnaf(i) =  - alphah_naf(k1) - betah_naf(k1)
       dfmkdr_dmkdr(i) =  - alpham_kdr(k1) - betam_kdr(k1)
       dfmka_dmka(i)  =   - alpham_ka (k1) - betam_ka (k1)
       dfhka_dhka(i)  =   - alphah_ka (k1) - betah_ka (k1)
       dfmk2_dmk2(i)  =   - alpham_k2 (k1) - betam_k2 (k1)
       dfhk2_dhk2(i)  =   - alphah_k2 (k1) - betah_k2 (k1)
       dfmkm_dmkm(i)  =   - alpham_km (k1) - betam_km (k1)
       dfmkc_dmkc(i)  =   - alpham_kc (k1) - betam_kc (k1)
       dfmcat_dmcat(i) =  - alpham_cat(k1) - betam_cat(k1)
       dfhcat_dhcat(i) =  - alphah_cat(k1) - betah_cat(k1)
       dfmcaL_dmcaL(i) =  - alpham_caL(k1) - betam_caL(k1)
       dfmar_dmar(i)  =   - alpham_ar (k1) - betam_ar (k1)

          end do

       dt2 = 0.5d0 * dt * dt

        do i = 1, numcomp
          v(i,L) = v(i,L) + dt * fv(i)
           do j = 1, numcomp
        v(i,L) = v(i,L) + dt2 * dfv_dv(i,j) * fv(j)
           end do
        v(i,L) = v(i,L) + dt2 * ( dfv_dchi(i) * fchi(i)
     X          + dfv_dmnaf(i) * fmnaf(i)
     X          + dfv_dhnaf(i) * fhnaf(i)
     X          + dfv_dmkdr(i) * fmkdr(i)
     X          + dfv_dmka(i)  * fmka(i)
     X          + dfv_dhka(i)  * fhka(i)
     X          + dfv_dmk2(i)  * fmk2(i)
     X          + dfv_dhk2(i)  * fhk2(i)
     X          + dfv_dmkm(i)  * fmkm(i)
     X          + dfv_dmkc(i)  * fmkc(i)
     X          + dfv_dmkahp(i)* fmkahp(i)
     X          + dfv_dmcat(i)  * fmcat(i)
     X          + dfv_dhcat(i) * fhcat(i)
     X          + dfv_dmcaL(i) * fmcaL(i)
     X          + dfv_dmar(i)  * fmar(i) )

        chi(i,L) = chi(i,L) + dt * fchi(i) + dt2 *
     X   (dfchi_dchi(i) * fchi(i) + dfchi_dv(i) * fv(i))
        mnaf(i,L) = mnaf(i,L) + dt * fmnaf(i) + dt2 *
     X   (dfmnaf_dmnaf(i) * fmnaf(i) + dfmnaf_dv(i)*fv(i))
        hnaf(i,L) = hnaf(i,L) + dt * fhnaf(i) + dt2 *
     X   (dfhnaf_dhnaf(i) * fhnaf(i) + dfhnaf_dv(i)*fv(i))
        mkdr(i,L) = mkdr(i,L) + dt * fmkdr(i) + dt2 *
     X   (dfmkdr_dmkdr(i) * fmkdr(i) + dfmkdr_dv(i)*fv(i))
        mka(i,L) =  mka(i,L) + dt * fmka(i) + dt2 *
     X   (dfmka_dmka(i) * fmka(i) + dfmka_dv(i) * fv(i))
        hka(i,L) =  hka(i,L) + dt * fhka(i) + dt2 *
     X   (dfhka_dhka(i) * fhka(i) + dfhka_dv(i) * fv(i))
        mk2(i,L) =  mk2(i,L) + dt * fmk2(i) + dt2 *
     X   (dfmk2_dmk2(i) * fmk2(i) + dfmk2_dv(i) * fv(i))
        hk2(i,L) =  hk2(i,L) + dt * fhk2(i) + dt2 *
     X   (dfhk2_dhk2(i) * fhk2(i) + dfhk2_dv(i) * fv(i))
        mkm(i,L) =  mkm(i,L) + dt * fmkm(i) + dt2 *
     X   (dfmkm_dmkm(i) * fmkm(i) + dfmkm_dv(i) * fv(i))
        mkc(i,L) =  mkc(i,L) + dt * fmkc(i) + dt2 *
     X   (dfmkc_dmkc(i) * fmkc(i) + dfmkc_dv(i) * fv(i))
        mkahp(i,L) = mkahp(i,L) + dt * fmkahp(i) + dt2 *
     X (dfmkahp_dmkahp(i)*fmkahp(i) + dfmkahp_dchi(i)*fchi(i))
        mcat(i,L) =  mcat(i,L) + dt * fmcat(i) + dt2 *
     X   (dfmcat_dmcat(i) * fmcat(i) + dfmcat_dv(i) * fv(i))
        hcat(i,L) =  hcat(i,L) + dt * fhcat(i) + dt2 *
     X   (dfhcat_dhcat(i) * fhcat(i) + dfhcat_dv(i) * fv(i))
        mcaL(i,L) =  mcaL(i,L) + dt * fmcaL(i) + dt2 *
     X   (dfmcaL_dmcaL(i) * fmcaL(i) + dfmcaL_dv(i) * fv(i))
        mar(i,L) =   mar(i,L) + dt * fmar(i) + dt2 *
     X   (dfmar_dmar(i) * fmar(i) + dfmar_dv(i) * fv(i))
         end do


              end do



2001          CONTINUE


1000    CONTINUE
        END



C  SETS UP TABLES FOR RATE FUNCTIONS
       SUBROUTINE SCORT_SETUP_suppyrRS
     X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
     X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
     X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
     X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
     X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
     X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
     X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
     X    alpham_km , betam_km , dalpham_km , dbetam_km ,
     X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
     X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
     X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
     X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
     X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)
      INTEGER I,J,K
      real*8 minf, hinf, taum, tauh, V, Z, shift_hnaf,
     X  shift_mkdr,
     X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
     X   dbetam_naf(0:640),
     X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
     X   dbetah_naf(0:640),
     X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
     X   dbetam_kdr(0:640),
     X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
     X   dbetam_ka(0:640),
     X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
     X   dbetah_ka(0:640),
     X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
     X   dbetam_k2(0:640),
     X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
     X   dbetah_k2(0:640),
     X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
     X   dbetam_km(0:640),
     X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
     X   dbetam_kc(0:640),
     X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
     X   dbetam_cat(0:640),
     X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
     X   dbetah_cat(0:640),
     X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
     X   dbetam_caL(0:640),
     X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
     X   dbetam_ar(0:640)
C FOR VOLTAGE, RANGE IS -120 TO +40 MV (absol.), 0.25 MV RESOLUTION


       DO 1, I = 0, 640
          V = dble(I)
          V = (V / 4.d0) - 120.d0

c gNa
           minf = 1.d0/(1.d0 + dexp((-V-38.d0)/10.d0))
           if (v.le.-30.d0) then
            taum = .025d0 + .14d0*dexp((v+30.d0)/10.d0)
           else
            taum = .02d0 + .145d0*dexp((-v-30.d0)/10.d0)
           endif
c from principal c. data, Martina & Jonas 1997, tau x 0.5
c Note that minf about the same for interneuron & princ. cell.
           alpham_naf(i) = minf / taum
           betam_naf(i) = 1.d0/taum - alpham_naf(i)

            shift_hnaf =  0.d0
        hinf = 1.d0/(1.d0 +
     x     dexp((v + shift_hnaf + 62.9d0)/10.7d0))
        tauh = 0.15d0 + 1.15d0/(1.d0+dexp((v+37.d0)/15.d0))
c from princ. cell data, Martina & Jonas 1997, tau x 0.5
            alphah_naf(i) = hinf / tauh
            betah_naf(i) = 1.d0/tauh - alphah_naf(i)

          shift_mkdr = 0.d0
c delayed rectifier, non-inactivating
       minf = 1.d0/(1.d0+dexp((-v-shift_mkdr-29.5d0)/10.0d0))
            if (v.le.-10.d0) then
             taum = .25d0 + 4.35d0*dexp((v+10.d0)/10.d0)
            else
             taum = .25d0 + 4.35d0*dexp((-v-10.d0)/10.d0)
            endif
              alpham_kdr(i) = minf / taum
              betam_kdr(i) = 1.d0 /taum - alpham_kdr(i)
c from Martina, Schultz et al., 1998. See espec. Table 1.

c A current: Huguenard & McCormick 1992, J Neurophysiol (TCR)
            minf = 1.d0/(1.d0 + dexp((-v-60.d0)/8.5d0))
            hinf = 1.d0/(1.d0 + dexp((v+78.d0)/6.d0))
        taum = .185d0 + .5d0/(dexp((v+35.8d0)/19.7d0) +
     x                            dexp((-v-79.7d0)/12.7d0))
        if (v.le.-63.d0) then
         tauh = .5d0/(dexp((v+46.d0)/5.d0) +
     x                  dexp((-v-238.d0)/37.5d0))
        else
         tauh = 9.5d0
        endif
           alpham_ka(i) = minf/taum
           betam_ka(i) = 1.d0 / taum - alpham_ka(i)
           alphah_ka(i) = hinf / tauh
           betah_ka(i) = 1.d0 / tauh - alphah_ka(i)

c h-current (anomalous rectifier), Huguenard & McCormick, 1992
           minf = 1.d0/(1.d0 + dexp((v+75.d0)/5.5d0))
           taum = 1.d0/(dexp(-14.6d0 -0.086d0*v) +
     x                   dexp(-1.87 + 0.07d0*v))
           alpham_ar(i) = minf / taum
           betam_ar(i) = 1.d0 / taum - alpham_ar(i)

c K2 K-current, McCormick & Huguenard
             minf = 1.d0/(1.d0 + dexp((-v-10.d0)/17.d0))
             hinf = 1.d0/(1.d0 + dexp((v+58.d0)/10.6d0))
            taum = 4.95d0 + 0.5d0/(dexp((v-81.d0)/25.6d0) +
     x                  dexp((-v-132.d0)/18.d0))
            tauh = 60.d0 + 0.5d0/(dexp((v-1.33d0)/200.d0) +
     x                  dexp((-v-130.d0)/7.1d0))
             alpham_k2(i) = minf / taum
             betam_k2(i) = 1.d0/taum - alpham_k2(i)
             alphah_k2(i) = hinf / tauh
             betah_k2(i) = 1.d0 / tauh - alphah_k2(i)

c voltage part of C-current, using 1994 kinetics, shift 60 mV
              if (v.le.-10.d0) then
       alpham_kc(i) = (2.d0/37.95d0)*dexp((v+50.d0)/11.d0 -
     x                                     (v+53.5)/27.d0)
       betam_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)-alpham_kc(i)
               else
       alpham_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)
       betam_kc(i) = 0.d0
               endif

c high-threshold gCa, from 1994, with 60 mV shift & no inactivn.
            alpham_cal(i) = 1.6d0/(1.d0+dexp(-.072d0*(v-5.d0)))
            betam_cal(i) = 0.1d0 * ((v+8.9d0)/5.d0) /
     x          (dexp((v+8.9d0)/5.d0) - 1.d0)

c M-current, from plast.f, with 60 mV shift
        alpham_km(i) = .02d0/(1.d0+dexp((-v-20.d0)/5.d0))
        betam_km(i) = .01d0 * dexp((-v-43.d0)/18.d0)

c T-current, from Destexhe, Neubig et al., 1998
         minf = 1.d0/(1.d0 + dexp((-v-56.d0)/6.2d0))
         hinf = 1.d0/(1.d0 + dexp((v+80.d0)/4.d0))
         taum = 0.204d0 + .333d0/(dexp((v+15.8d0)/18.2d0) +
     x                  dexp((-v-131.d0)/16.7d0))
          if (v.le.-81.d0) then
         tauh = 0.333 * dexp((v+466.d0)/66.6d0)
          else
         tauh = 9.32d0 + 0.333d0*dexp((-v-21.d0)/10.5d0)
          endif
              alpham_cat(i) = minf / taum
              betam_cat(i) = 1.d0/taum - alpham_cat(i)
              alphah_cat(i) = hinf / tauh
              betah_cat(i) = 1.d0 / tauh - alphah_cat(i)

1        CONTINUE

         do  i = 0, 639

      dalpham_naf(i) = (alpham_naf(i+1)-alpham_naf(i))/.25d0
      dbetam_naf(i) = (betam_naf(i+1)-betam_naf(i))/.25d0
      dalphah_naf(i) = (alphah_naf(i+1)-alphah_naf(i))/.25d0
      dbetah_naf(i) = (betah_naf(i+1)-betah_naf(i))/.25d0
      dalpham_kdr(i) = (alpham_kdr(i+1)-alpham_kdr(i))/.25d0
      dbetam_kdr(i) = (betam_kdr(i+1)-betam_kdr(i))/.25d0
      dalpham_ka(i) = (alpham_ka(i+1)-alpham_ka(i))/.25d0
      dbetam_ka(i) = (betam_ka(i+1)-betam_ka(i))/.25d0
      dalphah_ka(i) = (alphah_ka(i+1)-alphah_ka(i))/.25d0
      dbetah_ka(i) = (betah_ka(i+1)-betah_ka(i))/.25d0
      dalpham_k2(i) = (alpham_k2(i+1)-alpham_k2(i))/.25d0
      dbetam_k2(i) = (betam_k2(i+1)-betam_k2(i))/.25d0
      dalphah_k2(i) = (alphah_k2(i+1)-alphah_k2(i))/.25d0
      dbetah_k2(i) = (betah_k2(i+1)-betah_k2(i))/.25d0
      dalpham_km(i) = (alpham_km(i+1)-alpham_km(i))/.25d0
      dbetam_km(i) = (betam_km(i+1)-betam_km(i))/.25d0
      dalpham_kc(i) = (alpham_kc(i+1)-alpham_kc(i))/.25d0
      dbetam_kc(i) = (betam_kc(i+1)-betam_kc(i))/.25d0
      dalpham_cat(i) = (alpham_cat(i+1)-alpham_cat(i))/.25d0
      dbetam_cat(i) = (betam_cat(i+1)-betam_cat(i))/.25d0
      dalphah_cat(i) = (alphah_cat(i+1)-alphah_cat(i))/.25d0
      dbetah_cat(i) = (betah_cat(i+1)-betah_cat(i))/.25d0
      dalpham_caL(i) = (alpham_cal(i+1)-alpham_cal(i))/.25d0
      dbetam_caL(i) = (betam_cal(i+1)-betam_cal(i))/.25d0
      dalpham_ar(i) = (alpham_ar(i+1)-alpham_ar(i))/.25d0
      dbetam_ar(i) = (betam_ar(i+1)-betam_ar(i))/.25d0
       end do
2      CONTINUE

         do i = 640, 640
      dalpham_naf(i) =  dalpham_naf(i-1)
      dbetam_naf(i) =  dbetam_naf(i-1)
      dalphah_naf(i) = dalphah_naf(i-1)
      dbetah_naf(i) = dbetah_naf(i-1)
      dalpham_kdr(i) =  dalpham_kdr(i-1)
      dbetam_kdr(i) =  dbetam_kdr(i-1)
      dalpham_ka(i) =  dalpham_ka(i-1)
      dbetam_ka(i) =  dbetam_ka(i-1)
      dalphah_ka(i) =  dalphah_ka(i-1)
      dbetah_ka(i) =  dbetah_ka(i-1)
      dalpham_k2(i) =  dalpham_k2(i-1)
      dbetam_k2(i) =  dbetam_k2(i-1)
      dalphah_k2(i) =  dalphah_k2(i-1)
      dbetah_k2(i) =  dbetah_k2(i-1)
      dalpham_km(i) =  dalpham_km(i-1)
      dbetam_km(i) =  dbetam_km(i-1)
      dalpham_kc(i) =  dalpham_kc(i-1)
      dbetam_kc(i) =  dbetam_kc(i-1)
      dalpham_cat(i) =  dalpham_cat(i-1)
      dbetam_cat(i) =  dbetam_cat(i-1)
      dalphah_cat(i) =  dalphah_cat(i-1)
      dbetah_cat(i) =  dbetah_cat(i-1)
      dalpham_caL(i) =  dalpham_caL(i-1)
      dbetam_caL(i) =  dbetam_caL(i-1)
      dalpham_ar(i) =  dalpham_ar(i-1)
      dbetam_ar(i) =  dbetam_ar(i-1)
       end do   

4000   END
c Don't use this one
c      SUBROUTINE DEEPBASK_SETUP
c    X   (alpham_naf, betam_naf, dalpham_naf, dbetam_naf,
c    X    alphah_naf, betah_naf, dalphah_naf, dbetah_naf,
c    X    alpham_kdr, betam_kdr, dalpham_kdr, dbetam_kdr,
c    X    alpham_ka , betam_ka , dalpham_ka , dbetam_ka ,
c    X    alphah_ka , betah_ka , dalphah_ka , dbetah_ka ,
c    X    alpham_k2 , betam_k2 , dalpham_k2 , dbetam_k2 ,
c    X    alphah_k2 , betah_k2 , dalphah_k2 , dbetah_k2 ,
c    X    alpham_km , betam_km , dalpham_km , dbetam_km ,
c    X    alpham_kc , betam_kc , dalpham_kc , dbetam_kc ,
c    X    alpham_cat, betam_cat, dalpham_cat, dbetam_cat,
c    X    alphah_cat, betah_cat, dalphah_cat, dbetah_cat,
c    X    alpham_caL, betam_caL, dalpham_caL, dbetam_caL,
c    X    alpham_ar , betam_ar , dalpham_ar , dbetam_ar)

c     INTEGER I,J,K

c     real*8 minf, hinf, taum, tauh, V, Z, shift_hnaf,
c    X  shift_mkdr,
c    X alpham_naf(0:640),betam_naf(0:640),dalpham_naf(0:640),
c    X   dbetam_naf(0:640),
c    X alphah_naf(0:640),betah_naf(0:640),dalphah_naf(0:640),
c    X   dbetah_naf(0:640),
c    X alpham_kdr(0:640),betam_kdr(0:640),dalpham_kdr(0:640),
c    X   dbetam_kdr(0:640),
c    X alpham_ka(0:640), betam_ka(0:640),dalpham_ka(0:640) ,
c    X   dbetam_ka(0:640),
c    X alphah_ka(0:640), betah_ka(0:640), dalphah_ka(0:640),
c    X   dbetah_ka(0:640),
c    X alpham_k2(0:640), betam_k2(0:640), dalpham_k2(0:640),
c    X   dbetam_k2(0:640),
c    X alphah_k2(0:640), betah_k2(0:640), dalphah_k2(0:640),
c    X   dbetah_k2(0:640),
c    X alpham_km(0:640), betam_km(0:640), dalpham_km(0:640),
c    X   dbetam_km(0:640),
c    X alpham_kc(0:640), betam_kc(0:640), dalpham_kc(0:640),
c    X   dbetam_kc(0:640),
c    X alpham_cat(0:640),betam_cat(0:640),dalpham_cat(0:640),
c    X   dbetam_cat(0:640),
c    X alphah_cat(0:640),betah_cat(0:640),dalphah_cat(0:640),
c    X   dbetah_cat(0:640),
c    X alpham_caL(0:640),betam_caL(0:640),dalpham_caL(0:640),
c    X   dbetam_caL(0:640),
c    X alpham_ar(0:640), betam_ar(0:640), dalpham_ar(0:640),
c    X   dbetam_ar(0:640)
C FOR VOLTAGE, RANGE IS -120 TO +40 MV (absol.), 0.25 MV RESOLUTION


c      DO 1, I = 0, 640
c         V = dble (I)
c         V = (V / 4.d0) - 120.d0

c gNa
c          minf = 1.d0/(1.d0 + dexp((-V-38.d0)/10.d0))
c          if (v.le.-30.d0) then
c           taum = .0125d0 + .1525d0*dexp((v+30.d0)/10.d0)
c          else
c           taum = .02d0 + .145d0*dexp((-v-30.d0)/10.d0)
c          endif
c from interneuron data, Martina & Jonas 1997, tau x 0.5
c          alpham_naf(i) = minf / taum
c          betam_naf(i) = 1.d0/taum - alpham_naf(i)

c           shift_hnaf =  0.d0
c       hinf = 1.d0/(1.d0 +
c    x     dexp((v + shift_hnaf + 58.3d0)/6.7d0))
c       tauh = 0.225d0 + 1.125d0/(1.d0+dexp((v+37.d0)/15.d0))
c from interneuron data, Martina & Jonas 1997, tau x 0.5
c           alphah_naf(i) = hinf / tauh
c           betah_naf(i) = 1.d0/tauh - alphah_naf(i)

c         shift_mkdr = 0.d0
c delayed rectifier, non-inactivating
c      minf = 1.d0/(1.d0+dexp((-v-shift_mkdr-27.d0)/11.5d0))
c           if (v.le.-10.d0) then
c            taum = .25d0 + 4.35d0*dexp((v+10.d0)/10.d0)
c           else
c            taum = .25d0 + 4.35d0*dexp((-v-10.d0)/10.d0)
c           endif
c             alpham_kdr(i) = minf / taum
c             betam_kdr(i) = 1.d0 /taum - alpham_kdr(i)
c from Martina, Schultz et al., 1998

c A current: Huguenard & McCormick 1992, J Neurophysiol (TCR)
c           minf = 1.d0/(1.d0 + dexp((-v-60.d0)/8.5d0))
c           hinf = 1.d0/(1.d0 + dexp((v+78.d0)/6.d0))
c       taum = .185d0 + .5d0/(dexp((v+35.8d0)/19.7d0) +
c    x                            dexp((-v-79.7d0)/12.7d0))
c       if (v.le.-63.d0) then
c        tauh = .5d0/(dexp((v+46.d0)/5.d0) +
c    x                  dexp((-v-238.d0)/37.5d0))
c       else
c        tauh = 9.5d0
c       endif
c          alpham_ka(i) = minf/taum
c          betam_ka(i) = 1.d0 / taum - alpham_ka(i)
c          alphah_ka(i) = hinf / tauh
c          betah_ka(i) = 1.d0 / tauh - alphah_ka(i)

c h-current (anomalous rectifier), Huguenard & McCormick, 1992
c          minf = 1.d0/(1.d0 + dexp((v+75.d0)/5.5d0))
c          taum = 1.d0/(dexp(-14.6d0 -0.086d0*v) +
c    x                   dexp(-1.87 + 0.07d0*v))
c          alpham_ar(i) = minf / taum
c          betam_ar(i) = 1.d0 / taum - alpham_ar(i)

c K2 K-current, McCormick & Huguenard
c            minf = 1.d0/(1.d0 + dexp((-v-10.d0)/17.d0))
c            hinf = 1.d0/(1.d0 + dexp((v+58.d0)/10.6d0))
c           taum = 4.95d0 + 0.5d0/(dexp((v-81.d0)/25.6d0) +
c    x                  dexp((-v-132.d0)/18.d0))
c           tauh = 60.d0 + 0.5d0/(dexp((v-1.33d0)/200.d0) +
c    x                  dexp((-v-130.d0)/7.1d0))
c            alpham_k2(i) = minf / taum
c            betam_k2(i) = 1.d0/taum - alpham_k2(i)
c            alphah_k2(i) = hinf / tauh
c            betah_k2(i) = 1.d0 / tauh - alphah_k2(i)

c voltage part of C-current, using 1994 kinetics, shift 60 mV
c             if (v.le.-10.d0) then
c      alpham_kc(i) = (2.d0/37.95d0)*dexp((v+50.d0)/11.d0 -
c    x                                     (v+53.5)/27.d0)
c      betam_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)-alpham_kc(i)
c              else
c      alpham_kc(i) = 2.d0*dexp((-v-53.5d0)/27.d0)
c      betam_kc(i) = 0.d0
c              endif
c Speed-up of C kinetics here.
c         alpham_kc(i) = 2.d0 * alpham_kc(i)
c          betam_kc(i) = 2.d0 *  betam_kc(i)

c high-threshold gCa, from 1994, with 60 mV shift & no inactivn.
c           alpham_cal(i) = 1.6d0/(1.d0+dexp(-.072d0*(v-5.d0)))
c           betam_cal(i) = 0.1d0 * ((v+8.9d0)/5.d0) /
c    x          (dexp((v+8.9d0)/5.d0) - 1.d0)

c M-current, from plast.f, with 60 mV shift
c       alpham_km(i) = .02d0/(1.d0+dexp((-v-20.d0)/5.d0))
c       betam_km(i) = .01d0 * dexp((-v-43.d0)/18.d0)

c T-current, from Destexhe et al., 1996, pg. 170
c        minf = 1.d0/(1.d0 + dexp((-v-52.d0)/7.4d0))
c        hinf = 1.d0/(1.d0 + dexp((v+80.d0)/5.d0))
c        taum = 1.d0 + .33d0/(dexp((v+27.d0)/10.d0) +
c    x                  dexp((-v-102.d0)/15.d0))
c        tauh = 28.3d0 +.33d0/(dexp((v+48.d0)/4.d0) +
c    x                     dexp((-v-407.d0)/50.d0))
c             alpham_cat(i) = minf / taum
c             betam_cat(i) = 1.d0/taum - alpham_cat(i)
c             alphah_cat(i) = hinf / tauh
c             betah_cat(i) = 1.d0 / tauh - alphah_cat(i)

c1        CONTINUE

c        do 2, i = 0, 639

c     dalpham_naf(i) = (alpham_naf(i+1)-alpham_naf(i))/.25d0
c     dbetam_naf(i) = (betam_naf(i+1)-betam_naf(i))/.25d0
c     dalphah_naf(i) = (alphah_naf(i+1)-alphah_naf(i))/.25d0
c     dbetah_naf(i) = (betah_naf(i+1)-betah_naf(i))/.25d0
c     dalpham_kdr(i) = (alpham_kdr(i+1)-alpham_kdr(i))/.25d0
c     dbetam_kdr(i) = (betam_kdr(i+1)-betam_kdr(i))/.25d0
c     dalpham_ka(i) = (alpham_ka(i+1)-alpham_ka(i))/.25d0
c     dbetam_ka(i) = (betam_ka(i+1)-betam_ka(i))/.25d0
c     dalphah_ka(i) = (alphah_ka(i+1)-alphah_ka(i))/.25d0
c     dbetah_ka(i) = (betah_ka(i+1)-betah_ka(i))/.25d0
c     dalpham_k2(i) = (alpham_k2(i+1)-alpham_k2(i))/.25d0
c     dbetam_k2(i) = (betam_k2(i+1)-betam_k2(i))/.25d0
c     dalphah_k2(i) = (alphah_k2(i+1)-alphah_k2(i))/.25d0
c     dbetah_k2(i) = (betah_k2(i+1)-betah_k2(i))/.25d0
c     dalpham_km(i) = (alpham_km(i+1)-alpham_km(i))/.25d0
c     dbetam_km(i) = (betam_km(i+1)-betam_km(i))/.25d0
c     dalpham_kc(i) = (alpham_kc(i+1)-alpham_kc(i))/.25d0
c     dbetam_kc(i) = (betam_kc(i+1)-betam_kc(i))/.25d0
c     dalpham_cat(i) = (alpham_cat(i+1)-alpham_cat(i))/.25d0
c     dbetam_cat(i) = (betam_cat(i+1)-betam_cat(i))/.25d0
c     dalphah_cat(i) = (alphah_cat(i+1)-alphah_cat(i))/.25d0
c     dbetah_cat(i) = (betah_cat(i+1)-betah_cat(i))/.25d0
c     dalpham_caL(i) = (alpham_cal(i+1)-alpham_cal(i))/.25d0
c     dbetam_caL(i) = (betam_cal(i+1)-betam_cal(i))/.25d0
c     dalpham_ar(i) = (alpham_ar(i+1)-alpham_ar(i))/.25d0
c     dbetam_ar(i) = (betam_ar(i+1)-betam_ar(i))/.25d0
c2      CONTINUE

c        do i = 640, 640
c     dalpham_naf(i) =  dalpham_naf(i-1)
c     dbetam_naf(i) =  dbetam_naf(i-1)
c     dalphah_naf(i) = dalphah_naf(i-1)
c     dbetah_naf(i) = dbetah_naf(i-1)
c     dalpham_kdr(i) =  dalpham_kdr(i-1)
c     dbetam_kdr(i) =  dbetam_kdr(i-1)
c     dalpham_ka(i) =  dalpham_ka(i-1)
c     dbetam_ka(i) =  dbetam_ka(i-1)
c     dalphah_ka(i) =  dalphah_ka(i-1)
c     dbetah_ka(i) =  dbetah_ka(i-1)
c     dalpham_k2(i) =  dalpham_k2(i-1)
c     dbetam_k2(i) =  dbetam_k2(i-1)
c     dalphah_k2(i) =  dalphah_k2(i-1)
c     dbetah_k2(i) =  dbetah_k2(i-1)
c     dalpham_km(i) =  dalpham_km(i-1)
c     dbetam_km(i) =  dbetam_km(i-1)
c     dalpham_kc(i) =  dalpham_kc(i-1)
c     dbetam_kc(i) =  dbetam_kc(i-1)
c     dalpham_cat(i) =  dalpham_cat(i-1)
c     dbetam_cat(i) =  dbetam_cat(i-1)
c     dalphah_cat(i) =  dalphah_cat(i-1)
c     dbetah_cat(i) =  dbetah_cat(i-1)
c     dalpham_caL(i) =  dalpham_caL(i-1)
c     dbetam_caL(i) =  dbetam_caL(i-1)
c     dalpham_ar(i) =  dalpham_ar(i-1)
c     dbetam_ar(i) =  dbetam_ar(i-1)
c      end do   

c      END

        SUBROUTINE multipolarMAJ
C BRANCHED ACTIVE DENDRITES
     X             (GL,GAM,GKDR,GKA,GKC,GKAHP,GK2,GKM,
     X              GCAT,GCAL,GNAF,GNAP,GAR,
     X    CAFOR,JACOB,C,BETCHI,NEIGH,NNUM)
c Conductances: leak gL, coupling g, delayed rectifier gKDR, A gKA,
c C gKC, AHP gKAHP, K2 gK2, M gKM, low thresh Ca gCAT, high thresh
c gCAL, fast Na gNAF, persistent Na gNAP, h or anom. rectif. gAR.
c Note VAR = equil. potential for anomalous rectifier.
c Soma = comp. 1; 4 dendrites each with 13 compartments, 6-comp. axon
c Drop "glc"-like terms, just using "gl"-like
c CAFOR corresponds to "phi" in Traub et al., 1994
c Consistent set of units: nF, mV, ms, nA, microS

        INTEGER, PARAMETER:: numcomp = 59
        REAL*8 C(numcomp),GL(numcomp),GAM(0:numcomp,0:numcomp)
        REAL*8 GNAF(numcomp),GCAT(numcomp)
        REAL*8 GKDR(numcomp),GKA(numcomp),GKC(numcomp)
        REAL*8 GKAHP(numcomp),GCAL(numcomp),GAR(numcomp)
        REAL*8 GK2(numcomp),GKM(numcomp),GNAP(numcomp)
        REAL*8 JACOB(numcomp,numcomp)
        REAL*8 RI_SD,RI_AXON,RM_SD,RM_AXON,CDENS
        INTEGER LEVEL(numcomp)
        REAL*8 GNAF_DENS(0:9), GCAT_DENS(0:9), GKDR_DENS(0:9)
        REAL*8 GKA_DENS(0:9), GKC_DENS(0:9), GKAHP_DENS(0:9)
        REAL*8 GCAL_DENS(0:9), GK2_DENS(0:9), GKM_DENS(0:9)
        REAL*8 GNAP_DENS(0:9), GAR_DENS(0:9)
        REAL*8 RES, RINPUTi, ELEN(numcomp)
        REAL*8 RSOMA, PI, BETCHI(numcomp), CAFOR(numcomp)
        REAL*8 RAD(numcomp), LEN(numcomp), GAM1, GAM2
        REAL*8 RIN, D(numcomp), AREA(numcomp), RI, Z
        INTEGER NEIGH(numcomp,5), NNUM(numcomp), i, j, k, it
C FOR ESTABLISHING TOPOLOGY OF COMPARTMENTS

        RI_SD = 250.d0
        RM_SD = 50000.d0
c       RM_SD = 25000.d0
        RI_AXON = 100.d0
        RM_AXON = 1000.d0
        CDENS = 0.9d0

        PI = 3.14159d0

        gnaf_dens(0) = 400.d0
        gnaf_dens(1) =  60.d0
        gnaf_dens(2) =  60.d0
        gnaf_dens(3) =  60.d0
        do i = 4, 9
c         gnaf_dens(i) = 60.d0
          gnaf_dens(i) = 30.d0
        end do

        gkdr_dens(0) = 400.d0
        gkdr_dens(1) = 100.d0
        gkdr_dens(2) = 100.d0
        gkdr_dens(3) = 100.d0
        do i = 4, 9
         gkdr_dens(i) = 30.d0
c        gkdr_dens(i) = 60.d0
        end do

        gnap_dens(0) = 0.d0
        do i = 1, 9
c         gnap_dens(i) = 0.005d0 * gnaf_dens(i)
          gnap_dens(i) = 0.500d0 * gnaf_dens(i)
        end do

        gcat_dens(0) = 0.d0
        do i = 1, 3
          gcat_dens(i) = 0.05d0
        end do
        do i = 4, 9
          gcat_dens(i) = 0.5d0
        end do

        gcal_dens(0) = 0.d0
        do i = 1, 3
          gcal_dens(i) = 0.5d0
c         gcal_dens(i) = 0.1d0
        end do
        do i = 4, 9
c         gcal_dens(i) = 0.5d0
          gcal_dens(i) = 2.5d0
        end do

        gka_dens(0) = 1.d0
        gka_dens(1) =  2.d0 ! 4 April 2021
        gka_dens(2) =  1.d0
        gka_dens(3) =  1.d0
        do i = 4, 9
         gka_dens(i) = 1.0d0
        end do

        gkc_dens(0) = 0.d0
        do i = 1, 9
c        gkc_dens(i) = 10.00d0
         gkc_dens(i) =  0.00d0
        end do
        gkc_dens(1) = 10.d0

        gkm_dens(0) = 8.d0
        do i = 1, 9
         gkm_dens(i) = 6.00d0
        end do

        gk2_dens(0) = .5d0
        do i = 1, 9
         gk2_dens(i) = 0.00d0
        end do

        gkahp_dens(0) = 0.d0
        do i = 1, 9
         gkahp_dens(i) = 0.120d0
        end do

        gar_dens(0) = 0.d0
        do i = 1, 9
         gar_dens(i) = 0.02d0
        end do

c       WRITE   (6,9988)
9988    FORMAT(2X,'I',4X,'NADENS',' CADENS(L)',' KDRDEN',' KAHPDE',
     X     ' KCDENS',' KADENS')
c       DO 9989, I = 0, 9
        DO I = 0, 9
c         WRITE (6,9990) I, gnaf_dens(i), gcaL_dens(i), gkdr_dens(i),
c    X  gkahp_dens(i), gkc_dens(i), gka_dens(i)
9990    FORMAT(2X,I2,2X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2)
        END DO
9989    CONTINUE


        level(1) = 1
        do i = 2, 41, 13
         level(i) = 2
        end do
        do i = 3, 42, 13
           level(i) = 3
           level(i+1) = 3
        end do
        do i = 5, 44, 13
           level(i) = 4
           level(i+1) = 4
           level(i+2) = 4
        end do
        do i = 8, 47, 13
           level(i) = 5
           level(i+1) = 5
           level(i+2) = 5
        end do
        do i = 11, 50, 13
           level(i) = 6
           level(i+1) = 7
           level(i+2) = 8
           level(i+3) = 9
        end do

        do i = 54, 59
         level(i) = 0
        end do

c connectivity of axon
        nnum(54) = 2
        nnum(55) = 3
        nnum(56) = 3
        nnum(58) = 3
        nnum(57) = 1
        nnum(59) = 1
         neigh(54,1) =  1
         neigh(54,2) = 55
         neigh(55,1) = 54
         neigh(55,2) = 56
         neigh(55,3) = 58
         neigh(56,1) = 55
         neigh(56,2) = 57
         neigh(56,3) = 58
         neigh(58,1) = 55
         neigh(58,2) = 56
         neigh(58,3) = 59
         neigh(57,1) = 56
         neigh(59,1) = 58

c connectivity of SD part
          nnum(1) = 5
          neigh(1,1) = 54
          neigh(1,2) =  2
          neigh(1,3) = 15
          neigh(1,4) = 28
          neigh(1,5) = 41

          do i = 2, 41, 13
           nnum(i) = 3
           neigh(i,1) = 1
           neigh(i,2) = i + 1
           neigh(i,3) = i + 2
          end do

          do i = 3, 42, 13
           nnum(i) = 4
           neigh(i,1) = i - 1
           neigh(i,2) = i + 1
           neigh(i,3) = i + 2
           neigh(i,4) = i + 3
          end do

          do i = 4, 43, 13
           nnum(i) = 3
           neigh(i,1) = i - 2
           neigh(i,2) = i - 1
           neigh(i,3) = i + 3
          end do

          do i = 5, 44, 13
           nnum(i) = 3
           neigh(i,1) = i - 2
           neigh(i,2) = i + 1
           neigh(i,3) = i + 3
          end do

          do i = 6, 45, 13
           nnum(i) = 3
            neigh(i,1) = i - 3
            neigh(i,2) = i - 1
            neigh(i,3) = i + 3
          end do

          do i = 7, 46, 13
           nnum(i) = 2
           neigh(i,1) = i - 3
           neigh(i,2) = i + 3
          end do

          do i = 8, 47, 13
           nnum(i) = 2
           neigh(i,1) = i - 3
           neigh(i,2) = i + 3
          end do

          do i = 9, 48, 13
           nnum(i) = 1
           neigh(i,1) = i - 3
          end do

          do i = 10, 49, 13
           nnum(i) = 1
           neigh(i,1) = i - 3
          end do

          do i = 11, 50, 13
           nnum(i) = 2
           neigh(i,1) = i - 3
           neigh(i,2) = i + 1
          end do

          do i = 12, 51, 13
           nnum(i) = 2
           neigh(i,1) = i - 1
           neigh(i,2) = i + 1
          end do

          do i = 13, 52, 13
           nnum(i) = 2
           neigh(i,1) = i - 1
           neigh(i,2) = i + 1
          end do

          do i = 14, 53, 13
           nnum(i) = 1
           neigh(i,1) = i - 1
          end do

c        DO 332, I = 1, 59
         DO I = 1, numcomp
c          WRITE(6,3330) I, NEIGH(I,1),NEIGH(I,2),NEIGH(I,3),NEIGH(I,4),
c    X NEIGH(I,5)
3330     FORMAT(2X,I5,I5,I5,I5,I5,I5)
         END DO
332      CONTINUE
c         DO 858, I = 1, 59
          DO I = 1, 59
c          DO 858, J = 1, NNUM(I)
           DO J = 1, NNUM(I)
            K = NEIGH(I,J)
            IT = 0
c           DO 859, L = 1, NNUM(K)
            DO L = 1, NNUM(K)
             IF (NEIGH(K,L).EQ.I) IT = 1
            END DO
859         CONTINUE
             IF (IT.EQ.0) THEN
c             WRITE(6,8591) I, K
8591          FORMAT(' ASYMMETRY IN NEIGH MATRIX ',I4,I4)
             ENDIF
           END DO
           END DO
858       CONTINUE

c length and radius of axonal compartments
          do i = 54, 59
            len(i) = 50.d0
          end do
c         rad(54) = 0.80d0
c         rad(55) = 0.7d0
          rad(54) = 0.70d0
          rad(55) = 0.6d0
          do i = 56, 59
           rad(i) = 0.5d0
          end do

c  length and radius of SD compartments
          len(1) = 20.d0
          rad(1) = 7.5d0

          do i = 2, 53
c          len(i) = 40.d0
           len(i) = 80.d0 ! 4 April 2021
          end do

          rad(2) =   1.06d0
          rad(3) =   rad(2) / 1.59d0
          rad(4) =   rad(2) / 1.59d0
          rad(5) =   rad(2) / 2.53d0
          rad(6) =   rad(2) / 2.53d0
          rad(7) =   rad(2) / 1.59d0
          rad(8) =   rad(2) / 2.53d0
          rad(9) =   rad(2) / 2.53d0
          rad(10) =  rad(2) / 1.59d0
          rad(11) =  rad(2) / 2.53d0
          rad(12) =  rad(2) / 2.53d0
          rad(13) =  rad(2) / 2.53d0
          rad(14) =  rad(2) / 2.53d0

          do i = 15, 53
           rad(i) = rad(i-13)
          end do

c       WRITE(6,919)
919     FORMAT('COMPART.',' LEVEL ',' RADIUS ',' LENGTH(MU)')
c       DO 920, I = 1, 59
c920      WRITE(6,921) I, LEVEL(I), RAD(I), LEN(I)
921     FORMAT(I3,5X,I2,3X,F6.2,1X,F6.1,2X,F4.3)

c       DO 120, I = 1, 59
        DO I = 1, numcomp
           if (level(i).le.1) then 
          AREA(I) = 2.d0 * PI * RAD(I) * LEN(I)
           else
          AREA(I) = 4.d0 * PI * RAD(I) * LEN(I)
           endif
C NO CORRECTION FOR CONTRIBUTION OF SPINES TO AREA
! - in original bask.f, but change that here
          K = LEVEL(I)
          C(I) = CDENS * AREA(I) * (1.D-8)

           if (k.ge.1) then
          GL(I) = (1.D-2) * AREA(I) / RM_SD
           else
          GL(I) = (1.D-2) * AREA(I) / RM_AXON
           endif

          GNAF(I) = GNAF_DENS(K) * AREA(I) * (1.D-5)
          GNAP(I) = GNAP_DENS(K) * AREA(I) * (1.D-5)
          GCAT(I) = GCAT_DENS(K) * AREA(I) * (1.D-5)
          GKDR(I) = GKDR_DENS(K) * AREA(I) * (1.D-5)
          GKA(I) = GKA_DENS(K) * AREA(I) * (1.D-5)
          GKC(I) = GKC_DENS(K) * AREA(I) * (1.D-5)
          GKAHP(I) = GKAHP_DENS(K) * AREA(I) * (1.D-5)
          GCAL(I) = GCAL_DENS(K) * AREA(I) * (1.D-5)
          GK2(I) = GK2_DENS(K) * AREA(I) * (1.D-5)
          GKM(I) = GKM_DENS(K) * AREA(I) * (1.D-5)
          GAR(I) = GAR_DENS(K) * AREA(I) * (1.D-5)
c above conductances should be in microS
         END DO
120           continue

         Z = 0.d0
c        DO 1019, I = 2, 53
         DO I = 2, 53
           Z = Z + AREA(I)
         END DO
1019     CONTINUE
c        WRITE(6,1020) Z
1020     FORMAT(2X,' TOTAL DENDRITIC AREA ',F7.0)

c       DO 140, I = 1, 59
        DO I = 1, numcomp
c       DO 140, K = 1, NNUM(I)
        DO K = 1, NNUM(I)
         J = NEIGH(I,K)
           if (level(i).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM1 =100.d0 * PI * RAD(I) * RAD(I) / ( RI * LEN(I) )

           if (level(j).eq.0) then
               RI = RI_AXON
           else
               RI = RI_SD
           endif
         GAM2 =100.d0 * PI * RAD(J) * RAD(J) / ( RI * LEN(J) )

         GAM(I,J) = 2.d0/( (1.d0/GAM1) + (1.d0/GAM2) )
         END DO
         END DO
140     CONTINUE
c gam computed in microS

c       DO 299, I = 1, 59
        DO I = 1, numcomp
299       BETCHI(I) = .05d0
        END DO
        BETCHI( 1) =  .02d0

c       DO 300, I = 1, 59
        DO I = 1, numcomp
c300     D(I) = 2.D-4
300     D(I) = 1.D-4
        END DO
c       DO 301, I = 1, 59
        DO I = 1, numcomp
c        IF (LEVEL(I).EQ.1) D(I) = 5.D-3
         IF (LEVEL(I).EQ.1) D(I) = 2.D-4
        END DO
301     CONTINUE
C  NOTE NOTE NOTE  (DIFFERENT FROM SWONG)


c      DO 160, I = 1, 59
       DO I = 1, numcomp
160     CAFOR(I) = 5200.d0 / (AREA(I) * D(I))
       END DO
C     NOTE CORRECTION

c       do 200, i = 1, 59
        do i = 1, numcomp
200     C(I) = 1000.d0 * C(I)
        end do
C     TO GO FROM MICROF TO NF.

c     DO 909, I = 1, 59
      DO I = 1, numcomp
       JACOB(I,I) = - GL(I)
c     DO 909, J = 1, NNUM(I)
      DO J = 1, NNUM(I)
         K = NEIGH(I,J)
         IF (I.EQ.K) THEN
c            WRITE(6,510) I
510          FORMAT(' UNEXPECTED SYMMETRY IN NEIGH ',I4)
         ENDIF
         JACOB(I,K) = GAM(I,K)
         JACOB(I,I) = JACOB(I,I) - GAM(I,K)
       END DO
       END DO
909   CONTINUE

c 15 Jan. 2001: make correction for c(i)
          do i = 1, numcomp
          do j = 1, numcomp
             jacob(i,j) = jacob(i,j) / c(i)
          end do
          end do

c      DO 500, I = 1, 59
       DO I = 1, numcomp
c       WRITE (6,501) I,C(I)
501     FORMAT(1X,I2,' C(I) = ',F7.4)
        END DO
500     CONTINUE
        END

