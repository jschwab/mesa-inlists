require 'mesa_script'

m1s = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
m2s = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85]

m1s.each do |m1|
  m2s.each do |m2|

    id = 'he_%4.2f_co_%4.2f' % [m1, m2]
    
    # create the inlist to load the final model
    Inlist.make_inlist('inlists/inlist_load_%s' % [id]) do
      set_initial_model_number true
      initial_model_number 0
      load_saved_model true
      saved_model_name 'heco_models/%s.mod' % [id]
    end
    
    print "%s inlist_load_%s inlist_z_0.006 inlist_empty\n" % [id, id]
    
  end
end


m1s = [0.30]
m2s = [0.60]
ss = [1.2e9, 1e9, 8e8]

ss.each do |s|
  m1s.each do |m1|
    m2s.each do |m2|

      id = 'he_%4.2f_co_%4.2f_s_%6.1e' % [m1, m2, s]

      # create the inlist to load the final model
      Inlist.make_inlist('inlists/inlist_load_%s' % [id]) do
        set_initial_model_number true
        initial_model_number 0
        load_saved_model true
        saved_model_name 'heco_models/%s.mod' % [id]
      end
      
      print "%s inlist_load_%s inlist_z_0.006 inlist_empty\n" % [id, id]

    end
  end
end

