require 'mesa_script'

def write_composition(mhe, mco)

  filename = 'compositions/he_%4.2f_co_%4.2f.dat' % [mhe, mco]
  
  File.open(filename, 'w') do |file|
    
    qcore = mco / (mco + mhe)
    
    eps = 0.005

    nzones = 4
    nisos = 10

    def line
      fmt = "%18.12e %18.12e %18.12e %18.12e %18.12e %18.12e %18.12e %18.12e %18.12e %18.12e %18.12e\n"
      fmt % [@q, @xh1, @xhe3, @xhe4, @xc12, @xn14, @xo16, @xo18, @xne20, @xne22, @xmg24]
    end
    
    file.write("%s %s\n" % [nzones, nisos])

    # envelope (A00, fC=1.2)
    @xh1 = 0
    @xhe3 = 0
    @xhe4 = 0.9371
    @xc12 = 0.0409
    @xn14 = 0.0150
    @xo16 = 0.0035
    @xo18 = 0.0035
    @xne20 = 0.00
    @xne22 = 0.00
    @xmg24 = 0.00

    @q = 0
    file.write(line())

    @q = 1-qcore-eps
    file.write(line())

    # core (50/50 C/O)
    @xh1 = 0
    @xhe3 = 0
    @xhe4 = 0
    @xc12 = 0.5
    @xn14 = 0
    @xo16 = 0.5
    @xo18 = 0
    @xne20 = 0
    @xne22 = 0
    @xmg24 = 0


    @q = 1-qcore+eps
    file.write(line())

    @q = 1
    file.write(line())

  end

  return filename
  
end

m1s = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
m2s = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85]

m1s.each do |m1|
  m2s.each do |m2|

    id = 'he_%4.2f_co_%4.2f' % [m1, m2]
    
    # create the pms model inlist
    Inlist.make_inlist('inlists/inlist_pms_%s' % [id]) do
      initial_mass m1 + m2
      save_model_filename 'models/pms_%s.mod' % [id]
    end

    # write the composition
    filename = write_composition(m1, m2)

    # create the composition inlist
    Inlist.make_inlist('inlists/inlist_composition_%s' % [id]) do
      saved_model_name 'models/pms_%s.mod' % [id]
      save_model_filename 'models/composition_%s.mod' % [id]
      relax_composition_filename filename
    end

    # create the entropy inlist
    Inlist.make_inlist('inlists/inlist_entropy_%s' % [id]) do
      saved_model_name 'models/composition_%s.mod' % [id]
      save_model_filename 'models/entropy_%s.mod' % [id]
      filename_for_profile_when_terminate 'models/entropy_%s.data' % [id]
      x_ctrl(1, m2)
    end

    # create the finalization inlist
    Inlist.make_inlist('inlists/inlist_finalize_%s' % [id]) do
      saved_model_name 'models/entropy_%s.mod' % [id]
      save_model_filename 'final_models/%s.mod' % [id]
    end

    # create the inlist to load the final model
    Inlist.make_inlist('inlists/inlist_load_%s' % [id]) do
      saved_model_name '../make-HeCO-models/final_models/%s.mod' % [id]
    end
    
    print id+"\n"
    
  end
end


m1s = [0.30]
m2s = [0.60]
ss = [1.2e9, 1e9, 8e8]

ss.each do |s|
  m1s.each do |m1|
    m2s.each do |m2|

      id = 'he_%4.2f_co_%4.2f_s_%6.1e' % [m1, m2, s]

      # create the pms model inlist
      Inlist.make_inlist('inlists/inlist_pms_%s' % [id]) do
        initial_mass m1 + m2
        save_model_filename 'models/pms_%s.mod' % [id]
      end

      # write the composition
      filename = write_composition(m1, m2)

      # create the composition inlist
      Inlist.make_inlist('inlists/inlist_composition_%s' % [id]) do
        saved_model_name 'models/pms_%s.mod' % [id]
        save_model_filename 'models/composition_%s.mod' % [id]
        relax_composition_filename filename
      end

      # create the entropy inlist
      Inlist.make_inlist('inlists/inlist_entropy_%s' % [id]) do
        saved_model_name 'models/composition_%s.mod' % [id]
        save_model_filename 'models/entropy_%s.mod' % [id]
        filename_for_profile_when_terminate 'models/entropy_%s.data' % [id]
        x_ctrl(1, m2)
        x_ctrl(2, s)
      end

      # create the finalization inlist
      Inlist.make_inlist('inlists/inlist_finalize_%s' % [id]) do
        saved_model_name 'models/entropy_%s.mod' % [id]
        save_model_filename 'final_models/%s.mod' % [id]
      end

      # create the inlist to load the final model
      Inlist.make_inlist('inlists/inlist_load_%s' % [id]) do
        saved_model_name '../make-HeCO-models/final_models/%s.mod' % [id]
      end

      print id+"\n"

    end
  end
end

