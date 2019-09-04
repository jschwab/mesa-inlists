require 'mesa_script'

zs = ['0.0006', '0.002', '0.006', '0.02']

zs.each do |z|

  Inlist.make_inlist('inlists/inlist_z_%s' % [z]) do

    initial_z z
    zams_filename 'zams/he_zams_z%s.data' % [z]
    Zbase z.to_f
    
  end

end
