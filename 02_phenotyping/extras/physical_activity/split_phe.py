import pandas as pd

all_phes = pd.read_table('ashleylab_phenotypes_activity_and_fitness.txt')

all_phes.columns = ['FID', 'IID', 'Overall_acceleration_average', 'Standard_deviation_acceleration_average', 'MondayAverage', 'TuesdayAverage', 'WednesdayAverage', 'ThursdayAverage', 'FridayAverage', 'SaturdayAverage', 'SundayAverage', 'Average_acceleration_00:00_00:59', 'Average_acceleration_01:00_01:59', 'Average_acceleration_02:00_02:59', 'Average_acceleration_03:00_03:59', 'Average_acceleration_04:00_04:59', 'Average_acceleration_05:00_05:59', 'Average_acceleration_06:00_06:59', 'Average_acceleration_07:00_07:59', 'Average_acceleration_08:00_08:59', 'Average_acceleration_09:00_09:59', 'Average_acceleration_10:00_10:59', 'Average_acceleration_11:00_11:59', 'Average_acceleration_12:00_12:59', 'Average_acceleration_13:00_13:59', 'Average_acceleration_14:00_14:59', 'Average_acceleration_15:00_15:59', 'Average_acceleration_16:00_16:59', 'Average_acceleration_17:00_17:59', 'Average_acceleration_18:00_18:59', 'Average_acceleration_19:00_19:59', 'Average_acceleration_20:00_20:59', 'Average_acceleration_21:00_21:59', 'Average_acceleration_22:00_22:59', 'Average_acceleration_23:00_23:59', 'Fraction_acceleration_<=_1_milli_gravities', 'Fraction_acceleration_<=_2_milli_gravities', 'Fraction_acceleration_<=_3_milli_gravities', 'Fraction_acceleration_<=_4_milli_gravities', 'Fraction_acceleration_<=_5_milli_gravities', 'Fraction_acceleration_<=_6_milli_gravities', 'Fraction_acceleration_<=_7_milli_gravities', 'Fraction_acceleration_<=_8_milli_gravities', 'Fraction_acceleration_<=_9_milli_gravities', 'Fraction_acceleration_<=_10_milli_gravities', 'Fraction_acceleration_<=_11_milli_gravities', 'Fraction_acceleration_<=_12_milli_gravities', 'Fraction_acceleration_<=_13_milli_gravities', 'Fraction_acceleration_<=_14_milli_gravities', 'Fraction_acceleration_<=_15_milli_gravities', 'Fraction_acceleration_<=_16_milli_gravities', 'Fraction_acceleration_<=_17_milli_gravities', 'Fraction_acceleration_<=_18_milli_gravities', 'Fraction_acceleration_<=_19_milli_gravities', 'Fraction_acceleration_<=_20_milli_gravities', 'Fraction_acceleration_<=_25_milli_gravities', 'Fraction_acceleration_<=_30_milli_gravities', 'Fraction_acceleration_<=_35_milli_gravities', 'Fraction_acceleration_<=_40_milli_gravities', 'Fraction_acceleration_<=_45_milli_gravities', 'Fraction_acceleration_<=_50_milli_gravities', 'Fraction_acceleration_<=_55_milli_gravities', 'Fraction_acceleration_<=_60_milli_gravities', 'Fraction_acceleration_<=_65_milli_gravities', 'Fraction_acceleration_<=_70_milli_gravities', 'Fraction_acceleration_<=_75_milli_gravities', 'Fraction_acceleration_<=_80_milli_gravities', 'Fraction_acceleration_<=_85_milli_gravities', 'Fraction_acceleration_<=_90_milli_gravities', 'Fraction_acceleration_<=_95_milli_gravities', 'Fraction_acceleration_<=_100_milli_gravities', 'Fraction_acceleration_<=_125_milli_gravities', 'Fraction_acceleration_<=_150_milli_gravities', 'Fraction_acceleration_<=_175_milli_gravities', 'Fraction_acceleration_<=_200_milli_gravities', 'Fraction_acceleration_<=_225_milli_gravities', 'Fraction_acceleration_<=_250_milli_gravities', 'Fraction_acceleration_<=_275_milli_gravities', 'Fraction_acceleration_<=_300_milli_gravities', 'Fraction_acceleration_<=_325_milli_gravities', 'Fraction_acceleration_<=_350_milli_gravities', 'Fraction_acceleration_<=_375_milli_gravities', 'Fraction_acceleration_<=_400_milli_gravities', 'Fraction_acceleration_<=_425_milli_gravities', 'Fraction_acceleration_<=_450_milli_gravities', 'Fraction_acceleration_<=_475_milli_gravities', 'Fraction_acceleration_<=_500_milli_gravities', 'Fraction_acceleration_<=_600_milli_gravities', 'Fraction_acceleration_<=_700_milli_gravities', 'Fraction_acceleration_<=_800_milli_gravities', 'Fraction_acceleration_<=_900_milli_gravities', 'Fraction_acceleration_<=_1000_milli_gravities', 'Fraction_acceleration_<=_1100_milli_gravities', 'Fraction_acceleration_<=_1200_milli_gravities', 'Fraction_acceleration_<=_1300_milli_gravities', 'Fraction_acceleration_<=_1400_milli_gravities', 'Fraction_acceleration_<=_1500_milli_gravities', 'Fraction_acceleration_<=_1600_milli_gravities', 'Fraction_acceleration_<=_1700_milli_gravities', 'Fraction_acceleration_<=_1800_milli_gravities', 'Fraction_acceleration_<=_1900_milli_gravities', 'Fraction_acceleration_<=_2000_milli_gravities', 'Duration_walks', 'Frequency_strenuous_sports_last_4_weeks', 'Frequency_of_walking_for_pleasure_in_last_4_weeks', 'Number_days_moderate_physical_activity', 'Duration_moderate_activity', 'Number_days_vigorous_physical_activity', 'Duration_vigorous_activity', 'Time_spent_outdoors_summer', 'Time_spent_outdoors_winter', 'Duration_walking_for_pleasure', 'Duration_strenuous_sports', 'Number_days_walked_more_than_10_minutes', 'Usual_walking_pace', 'Number_transition_states_10mg', 'Number_transition_states_25mg', 'Discrete_wavelet_transform_signal_magnitude_vector', 'Discrete_wavelet_transform_signal_magnitude_vector_1', 'RHR_x', 'RHR_y', 'RHR_x.1', 'RHR_y.1', 'JobInvolesMainlyWalkingOrStanding', 'Job_involves_manual_or_physical_work', 'Job_involves_shift_work']

rel_cols = [col for col in all_phes.columns if not (col in ['FID', 'IID'])]

with open('/oak/stanford/groups/mrivas/users/guhan/repos/wiki/ukbb/icdinfo/icdinfo.txt', 'r') as icd:
    icdinfo = {line.split()[2]:line.split()[0] for line in icd}

for col in rel_cols:
    gbe_id = icdinfo[col] if col in icdinfo else col
    all_phes[['FID', 'IID', col]].to_csv(gbe_id + '.phe', sep='\t', index=False, header=False)