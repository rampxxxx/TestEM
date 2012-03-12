#------------------------------------------------------------------------------#
#-------------------------- TEST_ME make include file--------------------------#
#------------------------------------------------------------------------------#
HOME  =/home/rampxxxx/pfg/TestEM/source
HPREC = $(HOME)/prcsn
HEMRC = $(HOME)/em_src
HFK   = $(HEMRC)/FK
HIZ   = $(HEMRC)/IZ
HLR2K = $(HEMRC)/lr2k
HGR   = $(HEMRC)/grandi
HGRHF = $(HEMRC)/grandiHF
HCRR  = $(HEMRC)/carro
HNYGR = $(HEMRC)/nygren
HPKFST= $(HEMRC)/PKFstewart
HRYMC = $(HEMRC)/rmcl
HTENT = $(HEMRC)/tentusscher
HGPV  = $(HEMRC)/gpv

LIB   =$(HOME)/lib
#------------ Fortran Intel ----------------------------------------------------
F95   = ifort
#COM = -g -warn all -check bounds -p -traceback
COM  = -O3 -xW -static
MODU= -module $(HOME)/mod
