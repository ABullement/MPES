# This .R file can be used to lift over the relevant sensitivity analysis reported in the ESA of the tutorial paper.
# Please note: This code alone will not function!

# Sensitivity analysis inputs
default_knots <- c(0, 3.39, 5.23, 6.58)
default_hr <- c(1, 0.1, 6, 35)
default_gp <- c(34, 30196, 23080)
s1 <- c(0, 2.48490665, 4.787491743, 6.396929655)
s2 <- c(0, 4.094344562, 5.480638923, 6.58)
s3 <- c(0, log(12), log(24), log(36))
s4 <- c(0, 5.23, 6.58)
s5 <- c(0, 3.39, 6.58)
s6 <- c(0, 3.39, 5.23, 6.12, 6.58)
s8 <- c(1, 0.1, 1, 35)
s9 <- c(1, 0.1, 3, 35)
s10 <- c(1, 0.1, 6, 11)
s11 <- c(1, 0.1, 6, 13)
s13 <- c(34, 302, 231)
s14 <- c(34, 3019600, 2308000)

# Sensitivity analysis run lines
mpes_output_s1 <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = default_gp,
                           hr_specs = c(1, 0.1, 6, 35), knot_loc = s1, iter_spec = 1000, gp_choice = TRUE)

mpes_output_s2 <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = default_gp,
                           hr_specs = c(1, 0.1, 6, 35), knot_loc = s2, iter_spec = 1000, gp_choice = TRUE)

mpes_output_s3 <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = c(34, 30196, 23080),
                           hr_specs = c(1, 0.1, 6, 35), knot_loc = s3, iter_spec = 1000, gp_choice = TRUE)

mpes_output_s4 <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = default_gp,
                           hr_specs = c(1, 0.1, 6, 35), knot_loc = s4, iter_spec = 1000, gp_choice = TRUE)

mpes_output_s5 <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = default_gp,
                           hr_specs = c(1, 0.1, 6, 35), knot_loc = s5, iter_spec = 1000, gp_choice = TRUE)

mpes_output_s6 <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = default_gp,
                           hr_specs = c(1, 0.1, 6, 35), knot_loc = s6, iter_spec = 1000, gp_choice = TRUE)

mpes_output_s7 <- run_mpes(hr_adjust = FALSE, gen_pop_adjust = TRUE, gp_specs = default_gp,
                           hr_specs = c(1, 0.1, 6, 10), knot_loc = default_knots, iter_spec = 1000, gp_choice = TRUE)

mpes_output_s8 <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = default_gp,
                           hr_specs = s8, knot_loc = default_knots, iter_spec = 1000, gp_choice = TRUE)

mpes_output_s9 <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = default_gp,
                           hr_specs = s9, knot_loc = default_knots, iter_spec = 1000, gp_choice = TRUE)

mpes_output_s10 <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = default_gp,
                            hr_specs = s10, knot_loc = default_knots, iter_spec = 1000, gp_choice = TRUE)

mpes_output_s11 <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = default_gp,
                            hr_specs = s11, knot_loc = default_knots, iter_spec = 1000, gp_choice = TRUE)

mpes_output_s12 <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = FALSE, gp_specs = default_gp,
                            hr_specs = c(1, 0.1, 6, 35), knot_loc = default_knots, iter_spec = 1000, gp_choice = TRUE)

mpes_output_s13 <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = s13,
                            hr_specs = c(1, 0.1, 6, 35), knot_loc = default_knots, iter_spec = 1000, gp_choice = TRUE)

mpes_output_s14 <- run_mpes(hr_adjust = TRUE, gen_pop_adjust = TRUE, gp_specs = s14,
                            hr_specs = c(1, 0.1, 6, 35), knot_loc = default_knots, iter_spec = 1000, gp_choice = TRUE)

mpes_output_s15 <- run_mpes(hr_adjust = FALSE, gen_pop_adjust = FALSE, gp_specs = default_gp,
                            hr_specs = c(1, 0.1, 6, 35), knot_loc = default_knots, iter_spec = 1000, gp_choice = TRUE)
