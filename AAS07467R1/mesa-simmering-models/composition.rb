require 'mesa_script'

x_string = "fiducial"
xc12 = 0.4
xne22 = 0.015
xna23 = 0.00015
xo16 = 1.0 - xc12 - xne22 - xna23

Inlist.make_inlist('inlists/inlist_composition_%s' % [x_string]) do

  save_model_when_terminate true
  save_model_filename "models/composition_%s.mod" % [x_string]

  accrete_given_mass_fractions true
  num_accretion_species 4
  
  accretion_species_id(1, 'c12')
  accretion_species_xa(1,  xc12)
  accretion_species_id(2, 'o16')
  accretion_species_xa(2,  xo16)
  accretion_species_id(3, 'ne22')
  accretion_species_xa(3,  xne22)
  accretion_species_id(4, 'na23')
  accretion_species_xa(4,  xna23)

  mass_change 0.0

end

# from MR16
x12 = 4.05e-1
x20 = 1.34e-3
x22 = 1.37e-2
x24 = 4.42e-4

# from P17
x21 = 3.74e-5
x23 = 1.42e-4
x25 = 3.84e-5
x27 = 5.60e-5


["toy", "toy-lowC"].each do |x_string|

  x12 = 0.405 if x_string == "toy"
  x12 = 0.300 if x_string == "toy-lowC"
  
  x16 = 1.0 - (x12 + x20 + x22 + x24) - (x21 + x23 + x25 + x27)
  p x16

  Inlist.make_inlist('inlists/inlist_composition_%s' % [x_string]) do

    save_model_when_terminate true
    save_model_filename "models/composition_%s.mod" % [x_string]

    accrete_given_mass_fractions true
    num_accretion_species 9

    accretion_species_id(1, 'c12')
    accretion_species_xa(1,  x12)
    accretion_species_id(2, 'o16')
    accretion_species_xa(2,  x16)
    accretion_species_id(3, 'ne22')
    accretion_species_xa(3,  x22)
    accretion_species_id(4, 'ne20')
    accretion_species_xa(4,  x20)
    accretion_species_id(5, 'ne21')
    accretion_species_xa(5,  x21)
    accretion_species_id(6, 'na23')
    accretion_species_xa(6,  x23)
    accretion_species_id(7, 'mg25')
    accretion_species_xa(7,  x25)
    accretion_species_id(8, 'al27')
    accretion_species_xa(8,  x27)
    accretion_species_id(9, 'mg24')
    accretion_species_xa(9,  x24)

    mass_change 0.0

  end

end
