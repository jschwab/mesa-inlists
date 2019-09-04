require 'mesa_script'

ms = [0.7, 0.75, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0, 1.05, 1.1, 1.2]
ms.each do |m|
  
  Inlist.make_inlist('inlists/inlist_mass_%.3f' % m) do

    initial_mass m
    
  end

end
