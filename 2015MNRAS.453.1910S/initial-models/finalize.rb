require_relative 'mesa_script'

ne_frac = 9.0 / 19.0

mdots = [1e-7, 3e-7, 1e-6, 3e-6, 1e-5]
mdots.each do |mdot|

  mdot_string = ("%3.0e" % mdot).gsub(/-0/, 'm')

  Inlist.make_inlist('inlists/inlist_finalize_%s' % [mdot_string]) do

    load_saved_model true
    saved_model_name "initial.mod"

    save_model_when_terminate true
    save_model_filename "final.mod"

    set_tau_factor true
    set_to_this_tau_factor 300

    change_net true
    new_net_name "ecapture.net"

    log_center_density_limit 9.4
    delta_lgRho_cntr_hard_limit 3e-3
    delta_lgRho_cntr_limit 1e-3

    net_rate_factor 0

    use_Type2_opacities true
    Zbase 0.02

    varcontrol_target 1e-3
    mass_change mdot
    accrete_same_as_surface true

  end
end
