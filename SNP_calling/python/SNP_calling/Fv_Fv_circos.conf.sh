
# MINIMUM CIRCOS CONFIGURATION







# Defines unit length for ideogram and tick spacing, referenced

# using "u" prefix, e.g. 10u

#chromosomes_units           = 1000000



# Show all chromosomes in karyotype file. By default, this is

# true. If you want to explicitly specify which chromosomes

# to draw, set this to 'no' and use the 'chromosomes' parameter.

# chromosomes_display_default = yes



# Chromosome name, size and color definition

karyotype = /home/connellj/fusarium_venenatum/circos_out/Fv_Fv/Fv_Fg_genome_edited.txt



<ideogram>

  <spacing>

    # spacing between ideograms is 0.5% of the image

    #default = 0.005r

    default = 0.001r

    <pairwise A3_5_contig_18 PH1_1>

      # spacing between contig1 and FoL_chromosome_1 is 4x of 0.5% (2%) of image

      # The angle of the ideogram is also edited in <image> below.

      spacing = 10r

    </pairwise>

    <pairwise PH1_2 PH1_1>

      spacing = 4r

    </pairwise>

    <pairwise PH1_3 PH1_2>

      #spacing = 0.005r

      spacing = 4r

    </pairwise>

    <pairwise PH1_4 PH1_3>

      spacing = 4r

    </pairwise>

    <pairwise A3_5_contig_85 PH1_4>

      spacing = 10r

    </pairwise>

  </spacing>



  # ideogram position, thickness and fill

  radius           = 0.90r

  thickness        = 30p

  fill             = yes



  stroke_thickness = 3

  stroke_color     = black



  # ideogram labels

  # <<include ideogram.label.conf>>

  show_label        = no



  # show labels only for contigs 1-16 and

  # use the chromosome name as the label, but replace "contig" with "FoC"

  label_format     = eval( var(idx) < 16? replace(var(chr),"contig_","FoC") : "")



  # 5% of inner radius outside outer ideogram radius

  label_radius = dims(ideogram,radius_outer) + 0.15r

  label_size        = 40

  label_font        = bold

  label_parallel    = yes





  # ideogram cytogenetic bands, if defined in the karyotype file

  # <<include bands.conf>>

</ideogram>



# image size, background color, angular position

# of first ideogram, transparency levels, output

# file and directory

#

# it is best to include these parameters from etc/image.conf

# and override any using param* syntax

#

# e.g.

#<image>

# <<include etc/image.conf>>

# radius* = 500

# </image>

<image>

  # override the default angle_offset of -90 defined in etc/image.conf

  angle_offset* = -90

  #radius* = 500

  <<include etc/image.conf>> # included from Circos distribution



</image>



# Specify which chromosomes will be drawn and their orientation

chromosomes_reverse = A3_5_contig_1, A3_5_contig_2, A3_5_contig_18, A3_5_contig_30, A3_5_contig_44, A3_5_contig_46, A3_5_contig_43, A3_5_contig_70, A3_5_contig_60, A3_5_contig_54, A3_5_contig_61, A3_5_contig_58, A3_5_contig_31, A3_5_contig_52, A3_5_contig_68, A3_5_contig_16, A3_5_contig_21, A3_5_contig_10, A3_5_contig_25, A3_5_contig_37, A3_5_contig_77, A3_5_contig_74, A3_5_contig_82, A3_5_contig_75, A3_5_contig_78, A3_5_contig_34, A3_5_contig_45, A3_5_contig_53, A3_5_contig_35, A3_5_contig_5, A3_5_contig_11, A3_5_contig_50, A3_5_contig_69, A3_5_contig_15, A3_5_contig_86, A3_5_contig_38, A3_5_contig_81, A3_5_contig_33, A3_5_contig_8, A3_5_contig_57, A3_5_contig_13, A3_5_contig_47, A3_5_contig_23, A3_5_contig_19, A3_5_contig_79, A3_5_contig_51, A3_5_contig_20, A3_5_contig_49, A3_5_contig_55, A3_5_contig_80, A3_5_contig_48, A3_5_contig_40, A3_5_contig_39, A3_5_contig_73, A3_5_contig_3, PH1_4, PH1_3, PH1_2, PH1_1

chromosomes_order = A3_5_contig_18, A3_5_contig_30, A3_5_contig_44, A3_5_contig_29, A3_5_contig_46, A3_5_contig_43, A3_5_contig_70, A3_5_contig_60, A3_5_contig_54, A3_5_contig_7, A3_5_contig_76, A3_5_contig_58, A3_5_contig_52, A3_5_contig_68, A3_5_contig_84, A3_5_contig_16, A3_5_contig_28, A3_5_contig_67, A3_5_contig_21, A3_5_contig_10, A3_5_contig_25, A3_5_contig_37, A3_5_contig_65, A3_5_contig_77, A3_5_contig_74, A3_5_contig_82, A3_5_contig_26, A3_5_contig_75, A3_5_contig_78, A3_5_contig_34, A3_5_contig_45, A3_5_contig_53, A3_5_contig_14, A3_5_contig_32, A3_5_contig_24, A3_5_contig_17, A3_5_contig_35, A3_5_contig_5, A3_5_contig_50, A3_5_contig_69, A3_5_contig_15, A3_5_contig_31, A3_5_contig_6, A3_5_contig_86, A3_5_contig_38, A3_5_contig_81, A3_5_contig_61, A3_5_contig_62, A3_5_contig_12, A3_5_contig_1, A3_5_contig_22, A3_5_contig_66, A3_5_contig_57, A3_5_contig_13, A3_5_contig_9, A3_5_contig_59, A3_5_contig_47, A3_5_contig_23, A3_5_contig_19, A3_5_contig_64, A3_5_contig_79, A3_5_contig_51, A3_5_contig_20, A3_5_contig_41, A3_5_contig_49, A3_5_contig_63, A3_5_contig_2, A3_5_contig_33, A3_5_contig_55, A3_5_contig_80, A3_5_contig_11, A3_5_contig_48, A3_5_contig_36, A3_5_contig_42, A3_5_contig_40, A3_5_contig_83, A3_5_contig_39, A3_5_contig_27, A3_5_contig_4, A3_5_contig_73, A3_5_contig_3, A3_5_contig_8, A3_5_contig_56, A3_5_contig_71, A3_5_contig_72, A3_5_contig_85, PH1_4, PH1_3, PH1_2, PH1_1

# RGB/HSV color definitions, color lists, location of fonts,

# fill patterns

<<include etc/colors_fonts_patterns.conf>> # included from Circos distribution



# debugging, I/O an dother system parameters

<<include etc/housekeeping.conf>> # included from Circos distribution



# Include ticks

<<include /home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/analysis/synteny/Fv_Fg_ticks.conf>>

# Include a 2D plot

<<include /home/armita/git_repos/emr_repos/scripts/fusarium_venenatum/analysis/synteny/Fv_Fg_2D_plot.conf>>