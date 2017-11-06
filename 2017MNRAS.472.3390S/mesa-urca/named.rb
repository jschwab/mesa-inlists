require 'mesa_script'

names = ['SQB15', 'SQB15+', 'T13', 'F15']

names.each do |name|

  Inlist.make_inlist('inlists/inlist_io_%s' % [name]) do

    load_saved_model true
    saved_model_name "models/final_%s.mod" % [name]
    
    log_directory 'LOGS-named'
    star_history_name '%s.data' % [name]
    history_interval 1
    terminal_interval 100
    write_profiles_flag false
    
    pgstar_flag false
    
  end

  Inlist.make_inlist('inlists/inlist_%s' % [name]) do

    mlt_option 'none'
    
  end
  
end


