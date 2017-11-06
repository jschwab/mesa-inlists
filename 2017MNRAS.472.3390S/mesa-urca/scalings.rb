require 'mesa_script'

logxs = [-3.00, -2.75, -2.50, -2.25, -2.00, -1.75, -1.50]


def do_one_x(x_string)

  Inlist.make_inlist('inlists-scalings/inlist_%s' % [x_string]) do

    load_saved_model true
    saved_model_name "models/final_%s.mod" % [x_string]

    log_directory 'LOGS-scalings'
    star_history_name '%s.data' % [x_string]
    history_interval 1
    write_profiles_flag false

    log_center_density_limit 9.4

  end

end

def do_one_mdot(mdot_string)

  mdot = mdot_string[0].to_f * 10.0 ** (-mdot_string[-1].to_f)
  
  Inlist.make_inlist('inlists-scalings/inlist_mdot_%s' % [mdot_string]) do

    load_saved_model true
    saved_model_name "models/final_logxna-2.00.mod"

    log_directory 'LOGS-scalings'
    star_history_name 'mdot-%s.data' % [mdot_string]
    history_interval 1
    write_profiles_flag false
    
    log_center_density_limit 9.4

    mass_change mdot

    pgstar_flag false

  end

end

['na', 'mg', 'both'].each do |id|
  logxs.each do |logx|
    do_one_x("logx%s%5.2f" % [id, logx])
  end
end

mdots = ['1em9', '3em9', '1em8', '3em8', '1em7', '3em7', '1em6', '3em6']
mdots.each do |mdot_s|
  do_one_mdot(mdot_s)
end

  
