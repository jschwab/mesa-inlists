require 'mesa_script'

def do_one()
  
  Inlist.make_inlist('inlist_io') do

    load_saved_model true
    saved_model_name "models/%s.mod" % [ARGV[0]]

    if ARGV[1] then
      id = "%s-%s" % ARGV[0..1]
    else
      id = ARGV[0]
    end

    print(id)

    log_directory 'LOGS-%s' % [id]
    history_interval 1

    save_model_when_terminate true
    save_model_filename 'models/%s-at-core-He-ign.mod' % [id]

    pgstar_flag true

  end

end

do_one()

  
