require 'mesa_script'

masses = [1.300,
          1.305,
          1.310,
          1.315,
          1.320,
          1.325,
          1.330,
          1.335,
          1.340,
          1.345,
          1.350,
          1.355,
          1.360,
          1.365,
          1.370,
          1.375,
          1.380,
          1.385,
          1.390,
          1.395,
          1.400,
          1.405,
          1.410,
          1.415,
          1.420,
          1.425,
          1.430,
          1.435,
          1.440,
          1.445,
          1.450,
          1.455,
          1.460,
          1.465,
          1.470,
          1.475,
          1.480,
          1.485,
          1.490,
          1.495,
          1.500,
         ]

compositions = ['ne', 'o', 'neo', 'si', 'sis', 'fe']

compositions.each do |composition|
  masses.each do |mass|

    mass_string = ('%5.3f' % mass).sub(/\./,'p')
    comp_string = composition
    
    # make the pre-main sequence inlist
    id = '%s_%s' % [mass_string, comp_string]
    Inlist.make_inlist('inlists/inlist_pms_%s' % id) do
      initial_mass mass
      save_model_filename 'models/pms_%s.mod' % id

      case composition
      when 'ne'
        num_accretion_species 1
        accretion_species_id 1, 'ne20'
        accretion_species_xa 1, 1.0
      when 'o'
        num_accretion_species 1
        accretion_species_id 1, 'o16'
        accretion_species_xa 1, 1.0
      when 'si'
        num_accretion_species 1
        accretion_species_id 1, 'si28'
        accretion_species_xa 1, 1.0
      when 'fe'
        num_accretion_species 1
        accretion_species_id 1, 'fe56'
        accretion_species_xa 1, 1.0
      when 'neo'
        num_accretion_species 2
        accretion_species_id 1, 'o16'
        accretion_species_xa 1, 0.5
        accretion_species_id 2, 'ne20'
        accretion_species_xa 2, 0.5
      when 'sis'
        num_accretion_species 2
        accretion_species_id 1, 'si28'
        accretion_species_xa 1, 0.5
        accretion_species_id 2, 's32'
        accretion_species_xa 2, 0.5
      else
        puts composition
      end
      
    end

    # make the later evolution inlist
    Inlist.make_inlist('inlists/inlist_evolve_%s' % id) do
      saved_model_name 'models/pms_%s.mod' % id
      log_directory 'output/%s' % id
      filename_for_profile_when_terminate 'output/%s/final.profile' % id
    end

    # make the output director
    begin
      Dir.mkdir('output/%s' % id)
    rescue Exception => e  
      puts e.message  
    end
    
  end
end
