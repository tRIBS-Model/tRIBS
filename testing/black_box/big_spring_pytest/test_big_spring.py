import numpy as np
import pandas as pd


def test_parallel_serial_match(setup_data):
    # this test is intended to evaluate if the serial and parallel results from big spring match

    print("\nTesting congruence between parallel and serial simulations")
    parallel_result, serial_result, _ = setup_data

    diff_mask = parallel_result != serial_result
    different_columns = diff_mask.any(axis=0)
    result = different_columns.index[different_columns].tolist()

    if result:
        flag = True  # List is not empty
    else:
        flag = False  # List is empty

    assert not flag, f"The following variables do not match between the parallel and serial simulation: {result}"


def test_parallel_waterbalance(setup_data):
    # This test is intended to calculate the difference between P - Q - ET - del Storage; For perfect conservation of
    # water/mass P - Q - ET - del Storage = 0. It's unlikely to reach exactly 0 due to numerical precision,
    # etc. so currently a threshold is set as an acceptable difference. Hopefully this difference can be decreased
    # with increasing performance.

    mrf, _, r_par = setup_data

    print("\nTesting Big Spring (parallel) Water Balance")

    wb_threshold = 10  # mm/year

    drainage_area = r_par.int_spatial_vars['VAr'].sum()
    weights = r_par.int_spatial_vars['VAr'].values / drainage_area
    porosity = np.sum(r_par.int_spatial_vars['ThetaS'].values * weights)
    years = float(r_par.options['runtime']['value']) / (365.25 * 24)

    # calculate changes in stores--note the assumption is that boundaries are closed except for streams
    # so sat and unsat fluxes = 0, also currently there is no average canopy storage variables so this portion of the
    # storage is ignored.

    unsat_mm = (mrf['MSMU'].iloc[-1] * mrf['MDGW'].iloc[-1] - mrf['MSMU'].iloc[0] *
                mrf['MDGW'].iloc[0]) * porosity
    sat_mm = (mrf['MDGW'].iloc[0] - mrf['MDGW'].iloc[-1]) * porosity
    can_swe_mm = 10 * (mrf['AvInSn'].iloc[-1] - mrf['AvInSn'].iloc[0])
    swe_mm = 10 * (mrf['AvSWE'].iloc[-1] - mrf['AvSWE'].iloc[0])
    p_mm = mrf['MAP'].sum()
    et_mm = mrf['MET'].sum() - 10 * (mrf['AvSnSub'].sum() + mrf['AvSnEvap'].sum() + mrf['AvInSu'].sum())
    q_mm = mrf['Srf'].sum() * 3600 * 1000 / drainage_area,

    del_s = unsat_mm + sat_mm + can_swe_mm + swe_mm
    loss = et_mm + q_mm
    test = (p_mm - del_s - loss) / years  # mm/yr

    # test avg. annual value
    print(f"Test metric (P - Loss - delS) = {test} mm/yr")

    assert abs(test) < wb_threshold


# def test_parallel_efficiency_match(setup_data):
#     # This test is set up to evaluate the model performance against observational data. Currently, it will likely fail
#     # as there are known issues with the observational data. Currently this test should be skipped.
#     mrf, _, r_par = setup_data
#
#     ref_kge = 1 - np.sqrt(2)
#
#     print(f"\nTesting Model Efficiency.\nThreshold = {ref_kge} from Koben et al. 2011")
#
#     # model results
#     qout = r_par.get_qout_results()
#     qout.set_index('Time', inplace=True)
#     qout_hour = qout.resample('h').mean()
#
#     # need to be updated so it's more portable
#     obs_data = '/Users/wr/Desktop/big_spring_pytest/data/obs/Discharge_MS-1.csv'
#
#     q_obs = pd.read_csv(obs_data, skiprows=24)
#
#     # convert date time
#     q_obs['Date-Time'] = pd.to_datetime(q_obs['Date-Time'], format='%m/%d/%y %H:%M')
#     q_obs['Date-Time'] = q_obs['Date-Time'] - pd.Timedelta(
#         days=1)  # somehow this data was offset by one day, verifed w/ precip data
#     q_obs.set_index('Date-Time', inplace=True)
#     # convert discharge from ft^3/s to m^3/s
#     conversion_factor = 0.0283168
#     q_obs['Q_cubic_m_sec'] = q_obs['Discharge'] * conversion_factor
#
#     # clean up date & remove questionable data
#     remove = [-5, -4, -3, 7, 11, 22]
#     q_obs_trim = q_obs[~q_obs['Grade'].isin(remove)]
#
#     # resample to hourly time steps
#     q_obs_hour = q_obs_trim.resample('h').mean()
#     q_obs_hour = q_obs_hour.loc[~np.isnan(q_obs_hour.Discharge)]
#
#     subset_qobs, subset_qout = q_obs_hour.align(qout_hour, join='inner', axis=0)
#
#     kge = r_par.kling_gupta_efficiency(subset_qout.Qstrm_m3s, subset_qobs.Q_cubic_m_sec)
#
#     print(f"Kling-Gupta efficiency = {kge}")
#
#     assert kge > ref_kge
