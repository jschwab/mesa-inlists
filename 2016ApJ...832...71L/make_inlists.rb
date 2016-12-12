require 'mesa_script'

# make lewis number inlists
ncrit_ratios = [0.10, 0.15, 0.20, 0.25]
ncrit_ratios.each do |nc|
  lewis_numbers = [0.1, 0.3, 0.5, 1.0, 3.0, 5.0, 10.0]
  lewis_numbers.each do |le|

    le_string = ('%4.1f' % le).strip
    nc_string = ('%5.3f' % nc).strip
    id = "Le-%s_Nc-%s" % [le_string, nc_string]
    Inlist.make_inlist('inlists/inlist_%s' % id) do
      steps_to_take_before_terminate 1000
      use_other_d_mix true
      x_ctrl 1, nc
      x_logical_ctrl 1, true
      x_ctrl 2, le

      Profile_Panels1_title 'Le=%s, Nc=' % [le_string, nc_string]
      History_Panels1_title 'Le=%s, Nc=' % [le_string, nc_string]

    end
  end
end

# make constant (for N < Nc) mixing inlists
ncrit_ratios = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5]
ncrit_ratios.each do |nc|

    le_string = ('%5.3f' % nc).strip
    
    id = le_string
    Inlist.make_inlist('inlists/inlist_ncrit_%s' % id) do
      steps_to_take_before_terminate 1000
      use_other_d_mix true
      x_ctrl 1, nc
      x_logical_ctrl 1, false
      x_ctrl 2, 1e12

      Profile_Panels1_title 'Ncrit: Nc=%s Np' % le_string
      History_Panels1_title 'Ncrit: Nc=%s Np' % le_string

    end
end
