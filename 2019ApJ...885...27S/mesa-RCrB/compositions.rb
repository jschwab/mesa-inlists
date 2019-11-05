require 'mesa_script'

#     Zi   element   (ni/nz)sun         (Xi/Z)sun         <Ai>
#      6     C      2.451814E-01     1.709837E-01     1.201114E+01
#      7     N      6.158679E-02     5.008508E-02     1.400676E+01
#      8     O      5.005962E-01     4.650221E-01     1.599938E+01

xzC = 1.709837E-01
xzN = 5.008508E-02
xzO = 4.650221E-01

fN = 1.7
fO = 0.4
fCO = 1.1

zs = ['0.0006', '0.002', '0.006', '0.02']
#fCOs = [0.4, 0.7, 1.1, 1.5]

zs = ['0.006']
fCOs = [0.4, 0.8, 1.2, 1.6]

zs.each do |z|

  znum = z.to_f
  
  fCOs.each do |fCO|

    fC = fCO + fO
    xc = xzC * znum * 10**fC
    xn = xzN * znum * 10**fN
    xo = xzO * znum * 10**fO

    # put everything else in helium
    y = 1 - xc - xn - xo

    xc_s = '%5.3f' % [xc]
    
    Inlist.make_inlist('inlists/inlist_z_%s_A00_XC_%s' % [z,xc_s]) do

      kappa_lowT_prefix 'AESOPUS_GS98_RCrB_Z0.006'
      
      num_accretion_species 5

      accretion_species_id 1, 'he4'
      accretion_species_xa 1, y

      accretion_species_id 2, 'c12'
      accretion_species_xa 2,  xc

      accretion_species_id 3, 'n14'
      accretion_species_xa 3,  xn

      accretion_species_id 4, 'o16'
      accretion_species_xa 4,  0.5*xo

      accretion_species_id 5, 'o18'
      accretion_species_xa 5,  0.5*xo
      
    end

  end
end
