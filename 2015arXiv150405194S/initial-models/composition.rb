require_relative 'mesa_script'

ne_frac = 9.0 / 19.0

xmgs = [0.01, 0.02, 0.03, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
xmgs.each do |xmg|

  # xne = (ne_frac * (1.0 - xmg)).round(2)
  # xo = (1.0 - xmg - xne).round(2)

  xne = xmg
  xmg = 0.05
  xo = (1.0 - xmg - xne).round(2)
  
  x_string = "%02i%02i%02i" % [xo * 100, xne * 100, xmg * 100]

  Inlist.make_inlist('inlists/inlist_composition_%s' % [x_string]) do

      load_saved_model true
      saved_model_name 'models/pms.mod'

      save_model_when_terminate true
      save_model_filename "models/composition_%s.mod" % [x_string]

      relax_initial_to_xaccrete true
      num_steps_to_relax_composition 100

      set_tau_factor true
      set_to_this_tau_factor 300

      change_net true
      new_net_name 'ecapture.net'

      log_center_density_limit 7
      net_rate_factor 0

      use_Type2_opacities true
      Zbase 0.02

      accrete_given_mass_fractions true
      num_accretion_species 3
      accretion_species_id(1, 'o16')
      accretion_species_xa(1,  xo)
      accretion_species_id(2, 'ne20')
      accretion_species_xa(2, xne)
      accretion_species_id(3, 'mg24')
      accretion_species_xa(3, xmg)

      mass_change 0.0

  end
end
