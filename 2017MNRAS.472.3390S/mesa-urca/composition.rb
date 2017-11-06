require 'mesa_script'

logxs = [-3.00, -2.75, -2.50, -2.25, -2.00, -1.75, -1.50]

logxnas = logxs
logxnas.each do |logxna|

  xna = 10**logxna
  xne = 0.5 - xna
  xo = 0.5
  
  x_string = "logxna%5.2f" % logxna

  Inlist.make_inlist('inlists/inlist_composition_%s' % [x_string]) do

      save_model_when_terminate true
      save_model_filename "models/composition_%s.mod" % [x_string]

      accrete_given_mass_fractions true
      num_accretion_species 3
      accretion_species_id(1, 'o16')
      accretion_species_xa(1,  xo)
      accretion_species_id(2, 'ne20')
      accretion_species_xa(2, xne)
      accretion_species_id(3, 'na23')
      accretion_species_xa(3, xna)

      mass_change 0.0

  end

end

logxmgs = logxs
logxmgs.each do |logxmg|

  xmg = 10**logxmg
  xne = 0.5 - xmg
  xo = 0.5
  
  x_string = "logxmg%5.2f" % logxmg

  Inlist.make_inlist('inlists/inlist_composition_%s' % [x_string]) do

      save_model_when_terminate true
      save_model_filename "models/composition_%s.mod" % [x_string]

      accrete_given_mass_fractions true
      num_accretion_species 3
      accretion_species_id(1, 'o16')
      accretion_species_xa(1,  xo)
      accretion_species_id(2, 'ne20')
      accretion_species_xa(2, xne)
      accretion_species_id(3, 'mg25')
      accretion_species_xa(3, xmg)

      mass_change 0.0

  end

end

logxboths = logxs
logxboths.each do |logxboth|

  xboth = 10**logxboth
  xne = 0.5 - 2 * xboth
  xo = 0.5
  
  x_string = "logxboth%5.2f" % logxboth

  Inlist.make_inlist('inlists/inlist_composition_%s' % [x_string]) do

      save_model_when_terminate true
      save_model_filename "models/composition_%s.mod" % [x_string]

      accrete_given_mass_fractions true
      num_accretion_species 4
      accretion_species_id(1, 'o16')
      accretion_species_xa(1,  xo)
      accretion_species_id(2, 'ne20')
      accretion_species_xa(2, xne)
      accretion_species_id(3, 'na23')
      accretion_species_xa(3, xboth)
      accretion_species_id(4, 'mg25')
      accretion_species_xa(4, xboth)

      mass_change 0.0

  end

end
