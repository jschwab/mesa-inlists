require 'mesa_script'

xna = 0.05
xmg24 = 0.05
xmg25 = 0.01

xcs = [0.0, 0.001, 0.003, 0.01, 0.02, 0.03, 0.04, 0.05, 0.10]

xcs.each do |xc|

  # do models without Urca isotopes

  xne = 0.5 - (xc + xmg24)
  xo = 0.5

  x_string = "nourca_xc-%5.3f" % xc

  Inlist.make_inlist('inlists/inlist_composition_%s' % [x_string]) do

      save_model_when_terminate true
      save_model_filename "models/composition_%s.mod" % [x_string]

      accrete_given_mass_fractions true
      num_accretion_species 4
      accretion_species_id(1, 'o16')
      accretion_species_xa(1,  xo)
      accretion_species_id(2, 'ne20')
      accretion_species_xa(2, xne)
      accretion_species_id(3, 'mg24')
      accretion_species_xa(3, xmg24)
      accretion_species_id(4, 'c12')
      accretion_species_xa(4, xc)

      mass_change 0.0

  end

  Inlist.make_inlist('inlists/inlist_io_composition_%s' % [x_string]) do

      load_saved_model true
      saved_model_name "models/final_%s.mod" % [x_string]

      log_directory 'LOGS-carbon'
      star_history_name '%s.data' % [x_string]

      history_interval 1
      write_profiles_flag false

      pgstar_flag false

  end

end


xcs.each do |xc|

  xne = 0.5 - (xna + xmg24 + xmg25 + xc)
  xo = 0.5


  # do models with Urca isotopes

  x_string = "withurca_xc-%5.3f" % xc

  Inlist.make_inlist('inlists/inlist_composition_%s' % [x_string]) do

      save_model_when_terminate true
      save_model_filename "models/composition_%s.mod" % [x_string]

      accrete_given_mass_fractions true
      num_accretion_species 6
      accretion_species_id(1, 'o16')
      accretion_species_xa(1,  xo)
      accretion_species_id(2, 'ne20')
      accretion_species_xa(2, xne)
      accretion_species_id(3, 'mg24')
      accretion_species_xa(3, xmg24)
      accretion_species_id(4, 'na23')
      accretion_species_xa(4, xna)
      accretion_species_id(5, 'mg25')
      accretion_species_xa(5, xmg25)
      accretion_species_id(6, 'c12')
      accretion_species_xa(6, xc)

      mass_change 0.0

  end

  Inlist.make_inlist('inlists/inlist_io_composition_%s' % [x_string]) do

      load_saved_model true
      saved_model_name "models/final_%s.mod" % [x_string]

      log_directory 'LOGS-carbon'
      star_history_name '%s.data' % [x_string]

      history_interval 1
      write_profiles_flag false

      pgstar_flag false

  end

end
  
