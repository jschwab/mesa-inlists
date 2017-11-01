require 'mesa_script'

def do_one_mdot(mdot_string)

  mdot = mdot_string[0].to_f * 10.0 ** (-mdot_string[-1].to_f)
  
  Inlist.make_inlist('inlists/inlist_mdot_%s' % [mdot_string]) do

    load_saved_model true
    saved_model_name "models/final_SQB15+_small_3em8.mod"

    save_model_when_terminate true
    save_model_filename 'LOGS-mdots/final-mdot-%s.mod' % [mdot_string]

    write_profile_when_terminate true
    filename_for_profile_when_terminate 'LOGS-mdots/final-mdot-%s.data' % [mdot_string]
    
    log_directory 'LOGS-mdots'
    star_history_name 'mdot-%s.data' % [mdot_string]
    history_interval 1
    profile_interval 100000
    
    mass_change mdot

    pgstar_flag true

  end

end

mdots = ['1em9', '3em9', '1em8', '3em8', '1em7', '3em7', '1em6', '3em6']
mdots.each do |mdot_s|
  do_one_mdot(mdot_s)
end

  
