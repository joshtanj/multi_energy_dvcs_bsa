package multi_energy_dvcs_bsa;

import java.util.ArrayList;

import javax.swing.JFrame;

import org.jlab.clas.physics.LorentzVector;
import org.jlab.geom.prim.Vector3D;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.utils.groups.IndexedList;

public class rgk_6535MeV_p_e_correction_exe {
	
	static HipoDataSource readerelastic = new HipoDataSource();
	static HipoDataSource readerpiplusn = new HipoDataSource();
	static HipoDataSource readerdvcs = new HipoDataSource();
	
	static H2F helastic_theta_vs_p_e = new H2F("elastic_theta_vs_p_e", "elastic_theta_vs_p_e", 200, 0, 8, 200, 0, 45);
	static H2F helastic_theta_vs_phi_e = new H2F("elastic_theta_vs_phi_e", "elastic_theta_vs_phi_e", 200, -180, 180, 200, 0, 45);
	static H1F hW = new H1F("hW", "hW", 100, 0.7, 1.2);
	static H1F hpcorW = new H1F("hpcorW", "hpcorW", 100, 0.7, 1.2);
	
	static H2F hpiplusn_theta_vs_p_e = new H2F("piplusn_theta_vs_p_e", "piplusn_theta_vs_p_e", 200, 0, 8, 200, 0, 45);
	static H2F hpiplusn_theta_vs_phi_e = new H2F("piplusn_theta_vs_phi_e", "piplusn_theta_vs_phi_e", 200, -180, 180, 200, 0, 45);
	static H1F hmissing_n = new H1F("hmissing_n", "hmissing_n", 100, 0.7, 1.2);
	static H1F hpcor_missing_n = new H1F("hpcormissing_n", "hpcormissing_n", 100, 0.7, 1.2);
	static H2F hpiplusn_theta_vs_p_e_psreg1 = new H2F("piplusn_theta_vs_p_e_psreg1", "piplusn_theta_vs_p_e_psreg1", 200, 0, 8, 200, 0, 45);
	static H2F hpiplusn_theta_vs_phi_e_psreg1 = new H2F("piplusn_theta_vs_phi_e_psreg1", "piplusn_theta_vs_phi_e_psreg1",
															200, -180, 180, 200, 0, 45);
	static H1F hmissing_n_psreg1 = new H1F("hmissing_n_psreg1", "hmissing_n_psreg1", 100, 0.7, 1.2);
	static H1F hpcor_missing_n_psreg1 = new H1F("hpcormissing_n_psreg1", "hpcormissing_n_psreg1", 100, 0.7, 1.2);
	static H2F hpiplusn_theta_vs_p_e_psreg2 = new H2F("piplusn_theta_vs_p_e_psreg2", "piplusn_theta_vs_p_e_psreg2", 200, 0, 8, 200, 0, 45);
	static H2F hpiplusn_theta_vs_phi_e_psreg2 = new H2F("piplusn_theta_vs_phi_e_psreg2", "piplusn_theta_vs_phi_e_psreg2",
															200, -180, 180, 200, 0, 45);
	static H1F hmissing_n_psreg2 = new H1F("hmissing_n_psreg2", "hmissing_n_psreg2", 100, 0.7, 1.2);
	static H1F hpcor_missing_n_psreg2 = new H1F("hpcormissing_n_psreg2", "hpcormissing_n_psreg2", 100, 0.7, 1.2);
	
	static H2F hdvcs_theta_vs_p_e = new H2F("dvcs_theta_vs_p_e", "dvcs_theta_vs_p_e", 200, 0, 8, 200, 0, 45);
	static H2F hdvcs_theta_vs_phi_e = new H2F("dvcs_theta_vs_phi_e", "dvcs_theta_vs_phi_e", 200, -180, 180, 200, 0, 45);
	static H1F htheta_cone_gamma = new H1F("theta_cone_gamma", "theta_cone_gamma", 100, 0, 7);
	static H1F hX_eg_m = new H1F("X_eg_m", "X_eg_m", 100, 0, 3);
	static H1F hX_epg_pt = new H1F("X_epg_pt", "X_epg_pt", 100, 0, 0.8);
	static H1F hX_epg_E = new H1F("X_epg_E", "X_epg_E", 100, -2.25, 2.25);
	static H1F hpcor_theta_cone_gamma = new H1F("pcortheta_cone_gamma", "pcortheta_cone_gamma", 100, 0, 7);
	static H1F hpcor_X_eg_m = new H1F("pcorX_eg_m", "pcorX_eg_m", 100, 0, 3);
	static H1F hpcor_X_epg_pt = new H1F("pcorX_epg_pt", "pcorX_epg_pt", 100, 0, 0.8);
	static H1F hpcor_X_epg_E = new H1F("pcorX_epg_E", "pcorX_epg_E", 100, -2.25, 2.25);
	
	static IndexedList<H1F> histGroups_W_sec = new IndexedList<H1F>(1);
	static IndexedList<H1F> histGroups_p_cor_W_sec = new IndexedList<H1F>(1);
	
	public static void histos_W_sec() {
		for(int seci = 0; seci < 6; seci++) {
			H1F hW_sec = new H1F("hW_sec", "hW_sec", 100, 0.7, 1.2);
			histGroups_W_sec.add(hW_sec, seci);
			H1F hp_cor_W_sec = new H1F("hp_cor_W_sec", "hp_cor_W_sec", 100, 0.7, 1.2);
			histGroups_p_cor_W_sec.add(hp_cor_W_sec, seci);
		}
	}
	
	static IndexedList<H1F> histGroups_W_theta_phi_bin = new IndexedList<H1F>(3);
	static IndexedList<H1F> histGroups_p_cor_W_theta_phi_bin = new IndexedList<H1F>(3);
	
	public static void histos_W_vs_phi_bin() {
		for(int seci = 0; seci < 6; seci++) {
			for(int thetabin = 0; thetabin < 9; thetabin++) {
				for(int phibin = 0; phibin < 4; phibin++) {
					H1F hW_vs_phi_bin = new H1F("hW_vs_phi_bin", "hW_vs_phi_bin", 100, 0.7, 1.2);
					histGroups_W_theta_phi_bin.add(hW_vs_phi_bin, seci, thetabin, phibin);
					H1F hp_cor_W_theta_phi_bin = new H1F("hp_cor_W_theta_phi_bin", "hp_cor_W_theta_phi_bin", 100, 0.7, 1.2);
					histGroups_p_cor_W_theta_phi_bin.add(hp_cor_W_theta_phi_bin, seci, thetabin, phibin);
				}
			}
		}
	}
	
	static IndexedList<H1F> histGroups_missing_n_sec = new IndexedList<H1F>(1);
	static IndexedList<H1F> histGroups_p_cor_missing_n_sec = new IndexedList<H1F>(1);
	
	public static void histos_missing_n_sec() {
		for(int seci = 0; seci < 6; seci++) {
			H1F hmissing_n_sec = new H1F("hmissing_n_sec", "hmissing_n_sec", 100, 0.7, 1.2);
			histGroups_missing_n_sec.add(hmissing_n_sec, seci);
			H1F hp_cor_missing_n_sec = new H1F("hp_cor_missing_n_sec", "hp_cor_missing_n_sec", 100, 0.7, 1.2);
			histGroups_p_cor_missing_n_sec.add(hp_cor_missing_n_sec, seci);
		}
	}
	
	static IndexedList<H1F> histGroups_missing_n_psreg1_sec = new IndexedList<H1F>(1);
	static IndexedList<H1F> histGroups_p_cor_missing_n_psreg1_sec = new IndexedList<H1F>(1);
	
	public static void histos_missing_n_psreg1_sec() {
		for(int seci = 0; seci < 6; seci++) {
			H1F hmissing_n_psreg1_sec = new H1F("hmissing_n_psreg1_sec", "hmissing_n_psreg1_sec", 100, 0.7, 1.2);
			histGroups_missing_n_psreg1_sec.add(hmissing_n_psreg1_sec, seci);
			H1F hp_cor_missing_n_psreg1_sec = new H1F("hp_cor_missing_n_psreg1_sec", "hp_cor_missing_n_psreg1_sec", 100, 0.7, 1.2);
			histGroups_p_cor_missing_n_psreg1_sec.add(hp_cor_missing_n_psreg1_sec, seci);
		}
	}
	
	static IndexedList<H1F> histGroups_missing_n_psreg2_sec = new IndexedList<H1F>(1);
	static IndexedList<H1F> histGroups_p_cor_missing_n_psreg2_sec = new IndexedList<H1F>(1);
	
	public static void histos_missing_n_psreg2_sec() {
		for(int seci = 0; seci < 6; seci++) {
			H1F hmissing_n_psreg2_sec = new H1F("hmissing_n_psreg2_sec", "hmissing_n_psreg2_sec", 100, 0.7, 1.2);
			histGroups_missing_n_psreg2_sec.add(hmissing_n_psreg2_sec, seci);
			H1F hp_cor_missing_n_psreg2_sec = new H1F("hp_cor_missing_n_psreg2_sec", "hp_cor_missing_n_psreg2_sec", 100, 0.7, 1.2);
			histGroups_p_cor_missing_n_psreg2_sec.add(hp_cor_missing_n_psreg2_sec, seci);
		}
	}
	
	static IndexedList<H1F> histGroups_missing_n_theta_phi_bin = new IndexedList<H1F>(3);
	static IndexedList<H1F> histGroups_p_cor_missing_n_theta_phi_bin = new IndexedList<H1F>(3);
	
	public static void histos_missing_n_vs_phi_bin() {
		for(int seci = 0; seci < 6; seci++) {
			for(int thetabin = 0; thetabin < 9; thetabin++) {
				for(int phibin = 0; phibin < 4; phibin++) {
					H1F hmissing_n_vs_phi_bin = new H1F("hmissing_n_vs_phi_bin", "hmissing_n_vs_phi_bin", 100, 0.7, 1.2);
					histGroups_missing_n_theta_phi_bin.add(hmissing_n_vs_phi_bin, seci, thetabin, phibin);
					H1F hp_cor_missing_n_theta_phi_bin = new H1F("hp_cor_missing_n_theta_phi_bin", "hp_cor_missing_n_theta_phi_bin",
																	100, 0.5, 1.5);
					histGroups_p_cor_missing_n_theta_phi_bin.add(hp_cor_missing_n_theta_phi_bin, seci, thetabin, phibin);
				}
			}
		}
	}
	
	static IndexedList<H1F> histGroups_missing_n_theta_phi_bin_psreg1 = new IndexedList<H1F>(3);
	static IndexedList<H1F> histGroups_p_cor_missing_n_theta_phi_bin_psreg1 = new IndexedList<H1F>(3);
	
	public static void histos_missing_n_vs_phi_bin_psreg1() {
		for(int seci = 0; seci < 6; seci++) {
			for(int thetabin = 0; thetabin < 9; thetabin++) {
				for(int phibin = 0; phibin < 4; phibin++) {
					H1F hmissing_n_vs_phi_bin_psreg1 = new H1F("hmissing_n_vs_phi_bin_psreg1", "hmissing_n_vs_phi_bin_psreg1", 100, 0.7, 1.2);
					histGroups_missing_n_theta_phi_bin_psreg1.add(hmissing_n_vs_phi_bin_psreg1, seci, thetabin, phibin);
					H1F hp_cor_missing_n_theta_phi_bin_psreg1 = new H1F("hp_cor_missing_n_theta_phi_bin_psreg1",
																			"hp_cor_missing_n_theta_phi_bin_psreg1", 100, 0.7, 1.2);
					histGroups_p_cor_missing_n_theta_phi_bin_psreg1.add(hp_cor_missing_n_theta_phi_bin_psreg1, seci, thetabin, phibin);
				}
			}
		}
	}
	
	static IndexedList<H1F> histGroups_missing_n_theta_phi_bin_psreg2 = new IndexedList<H1F>(3);
	static IndexedList<H1F> histGroups_p_cor_missing_n_theta_phi_bin_psreg2 = new IndexedList<H1F>(3);
	
	public static void histos_missing_n_vs_phi_bin_psreg2() {
		for(int seci = 0; seci < 6; seci++) {
			for(int thetabin = 0; thetabin < 9; thetabin++) {
				for(int phibin = 0; phibin < 4; phibin++) {
					H1F hmissing_n_vs_phi_bin_psreg2 = new H1F("hmissing_n_vs_phi_bin_psreg2", "hmissing_n_vs_phi_bin_psreg2", 100, 0.7, 1.2);
					histGroups_missing_n_theta_phi_bin_psreg2.add(hmissing_n_vs_phi_bin_psreg2, seci, thetabin, phibin);
					H1F hp_cor_missing_n_theta_phi_bin_psreg2 = new H1F("hp_cor_missing_n_theta_phi_bin_psreg2",
																			"hp_cor_missing_n_theta_phi_bin_psreg2", 100, 0.7, 1.2);
					histGroups_p_cor_missing_n_theta_phi_bin_psreg2.add(hp_cor_missing_n_theta_phi_bin_psreg2, seci, thetabin, phibin);
				}
			}
		}
	}
	
	static IndexedList<H1F> histGroups_theta_cone_gamma_sec = new IndexedList<H1F>(1);
	static IndexedList<H1F> histGroups_p_cor_theta_cone_gamma_sec = new IndexedList<H1F>(1);
	static IndexedList<H1F> histGroups_X_eg_m_sec = new IndexedList<H1F>(1);
	static IndexedList<H1F> histGroups_p_cor_X_eg_m_sec = new IndexedList<H1F>(1);
	static IndexedList<H1F> histGroups_X_epg_pt_sec = new IndexedList<H1F>(1);
	static IndexedList<H1F> histGroups_p_cor_X_epg_pt_sec = new IndexedList<H1F>(1);
	static IndexedList<H1F> histGroups_X_epg_E_sec = new IndexedList<H1F>(1);
	static IndexedList<H1F> histGroups_p_cor_X_epg_E_sec = new IndexedList<H1F>(1);
	
	public static void histos_dvcs_ev_sec() {
		for(int seci = 0; seci < 6; seci++) {
			H1F htheta_cone_gamma_sec = new H1F("htheta_cone_gamma_sec", "htheta_cone_gamma_sec", 100, 0, 7 );
			histGroups_theta_cone_gamma_sec.add(htheta_cone_gamma_sec, seci);
			H1F hp_cor_theta_cone_gamma_sec = new H1F("hp_cor_theta_cone_gamma_sec", "hp_cor_theta_cone_gamma_sec", 100, 0, 7 );
			histGroups_p_cor_theta_cone_gamma_sec.add(hp_cor_theta_cone_gamma_sec, seci);
			H1F hX_eg_m_sec = new H1F("hX_eg_m_sec", "hX_eg_m_sec", 100, 0, 3);
			histGroups_X_eg_m_sec.add(hX_eg_m_sec, seci);
			H1F hp_cor_X_eg_m_sec = new H1F("hp_cor_X_eg_m_sec", "hp_cor_X_eg_m_sec", 100, 0, 3);
			histGroups_p_cor_X_eg_m_sec.add(hp_cor_X_eg_m_sec, seci);
			H1F hX_epg_pt_sec = new H1F("hX_epg_pt_sec", "hX_epg_pt_sec", 100, 0, 0.8);
			histGroups_X_epg_pt_sec.add(hX_epg_pt_sec, seci);
			H1F hp_cor_X_epg_pt_sec = new H1F("hp_cor_X_epg_pt_sec", "hp_cor_X_epg_pt_sec", 100, 0, 0.8);
			histGroups_p_cor_X_epg_pt_sec.add(hp_cor_X_epg_pt_sec, seci);
			H1F hX_epg_E_sec = new H1F("hX_epg_E_sec", "hX_epg_E_sec", 100, -2.25, 2.25);
			histGroups_X_epg_E_sec.add(hX_epg_E_sec, seci);
			H1F hp_cor_X_epg_E_sec = new H1F("hp_cor_X_epg_E_sec", "hp_cor_X_epg_E_sec", 100, -2.25, 2.25);
			histGroups_p_cor_X_epg_E_sec.add(hp_cor_X_epg_E_sec, seci);
		}
	}
	
	static IndexedList<H1F> histGroups_theta_cone_gamma_theta_phi_bin = new IndexedList<H1F>(3);
	static IndexedList<H1F> histGroups_p_cor_theta_cone_gamma_theta_phi_bin = new IndexedList<H1F>(3);
	static IndexedList<H1F> histGroups_X_eg_m_theta_phi_bin = new IndexedList<H1F>(3);
	static IndexedList<H1F> histGroups_p_cor_X_eg_m_theta_phi_bin = new IndexedList<H1F>(3);
	static IndexedList<H1F> histGroups_X_epg_pt_theta_phi_bin = new IndexedList<H1F>(3);
	static IndexedList<H1F> histGroups_p_cor_X_epg_pt_theta_phi_bin = new IndexedList<H1F>(3);
	static IndexedList<H1F> histGroups_X_epg_E_theta_phi_bin = new IndexedList<H1F>(3);
	static IndexedList<H1F> histGroups_p_cor_X_epg_E_theta_phi_bin = new IndexedList<H1F>(3);
	
	public static void histos_dvcs_ev_vs_phi_bin() {
		for(int seci = 0; seci < 6; seci++) {
			for(int thetabin = 0; thetabin < 9; thetabin++) {
				for(int phibin = 0; phibin < 4; phibin++) {
					H1F htheta_cone_gamma_vs_phi_bin = new H1F("htheta_cone_gamma_vs_phi_bin", "htheta_cone_gamma_vs_phi_bin", 100, 0, 7 );
					histGroups_theta_cone_gamma_theta_phi_bin.add(htheta_cone_gamma_vs_phi_bin, seci, thetabin, phibin);
					H1F hp_cor_theta_cone_gamma_theta_phi_bin = new H1F("hp_cor_theta_cone_gamma_theta_phi_bin",
																			"hp_cor_theta_cone_gamma_theta_phi_bin", 100, 0, 7 );
					histGroups_p_cor_theta_cone_gamma_theta_phi_bin.add(hp_cor_theta_cone_gamma_theta_phi_bin, seci, thetabin, phibin);
					H1F hX_eg_m_vs_phi_bin = new H1F("hX_eg_m_vs_phi_bin", "hX_eg_m_vs_phi_bin", 100, 0, 3);
					histGroups_X_eg_m_theta_phi_bin.add(hX_eg_m_vs_phi_bin, seci, thetabin, phibin);
					H1F hp_cor_X_eg_m_theta_phi_bin = new H1F("hp_cor_X_eg_m_theta_phi_bin", "hp_cor_X_eg_m_theta_phi_bin", 100, 0, 3);
					histGroups_p_cor_X_eg_m_theta_phi_bin.add(hp_cor_X_eg_m_theta_phi_bin, seci, thetabin, phibin);
					H1F hX_epg_pt_vs_phi_bin = new H1F("hX_epg_pt_vs_phi_bin", "hX_epg_pt_vs_phi_bin", 100, 0, 0.8);
					histGroups_X_epg_pt_theta_phi_bin.add(hX_epg_pt_vs_phi_bin, seci, thetabin, phibin);
					H1F hp_cor_X_epg_pt_theta_phi_bin = new H1F("hp_cor_X_epg_pt_theta_phi_bin", "hp_cor_X_epg_pt_theta_phi_bin", 100, 0, 0.8);
					histGroups_p_cor_X_epg_pt_theta_phi_bin.add(hp_cor_X_epg_pt_theta_phi_bin, seci, thetabin, phibin);
					H1F hX_epg_vs_phi_bin = new H1F("hX_epg_vs_phi_bin", "hX_epg_vs_phi_bin", 100, -2.25, 2.25);
					histGroups_X_epg_E_theta_phi_bin.add(hX_epg_vs_phi_bin, seci, thetabin, phibin);
					H1F hp_cor_X_epg_E_theta_phi_bin = new H1F("hp_cor_X_epg_E_theta_phi_bin", "hp_cor_X_epg_E_theta_phi_bin", 100, -2.25, 2.25);
					histGroups_p_cor_X_epg_E_theta_phi_bin.add(hp_cor_X_epg_E_theta_phi_bin, seci, thetabin, phibin);
				}
			}
		}
	}
	
	static double rc = Math.PI/180;

	static void processEventelastic(DataEvent event, int corEventCounter, IndexedList<double[]> thetaCorFac) {
		
		if(event.hasBank("REC::Particle") && event.hasBank("REC::Track") && event.hasBank("REC::Traj"))
		{
			HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
			int ecount = 0;
			int protoncount = 0;
			int e_pindex = -1;
			int proton_pindex = -1;
			float M = (float) 0.938;
			float E = (float) 6.535;
			float W = -20;
			float W_p_cor = -20;
			float theta_deg_e = 0;
			float phi_deg_e = -200;
			float phi_rot_deg_e = -200;
			double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0);
			LorentzVector p_proton = new LorentzVector(0, 0, 0, 0);
			Vector3D v_e = new Vector3D (0, 0, 0);
			Vector3D v_proton = new Vector3D (0, 0, 0);
			float cdchi2_proton = 100;
			short cdNDF_proton = 100;
			Vector3D dc_hit_e = new Vector3D (0, 0, 0);
			Vector3D dc_hit_rot_e = new Vector3D (0, 0, 0);
			double pcf = 0;
			for(int i = 0; i < rec.rows(); i++)
			{
				int pid = rec.getInt("pid", i);
				float px = rec.getFloat("px", i);
				float py = rec.getFloat("py", i);
				float pz = rec.getFloat("pz", i);
				float vx = rec.getFloat("vx", i);
				float vy = rec.getFloat("vy", i);
				float vz = rec.getFloat("vz", i);
				float p = (float) Math.sqrt((px*px)+(py*py)+(pz*pz));
				if(pid == 11)
				{
					ecount++;
					if(ecount == 1)
					{
						p_e = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.000511*0.000511)));
						v_e = new Vector3D (vx, vy, vz);
						e_pindex = i;
					}
				}
				if(pid == 2212)
				{
					protoncount++;
					if(protoncount == 1)
					{
						p_proton = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.938272*0.938272)));
						v_proton = new Vector3D (vx, vy, vz);
						proton_pindex = i;
					}
				}
			}
			HipoDataBank rectrac = (HipoDataBank) event.getBank("REC::Track");
			byte e_detector = 0;
			byte e_sector = 0;
			byte proton_detector = 0;
			for(int j = 0; j < rectrac.rows(); j++)
			{
				short pindex = rectrac.getShort("pindex", j);
				byte detector = rectrac.getByte("detector", j);
				byte sector = rectrac.getByte("sector", j);
				float chi2 = rectrac.getFloat("chi2", j);
				short NDF = rectrac.getShort("NDF", j);
				if(pindex == e_pindex)
				{
					e_detector = detector;
					e_sector = sector;
				}
				if(pindex == proton_pindex)
				{
					proton_detector = detector;
					cdchi2_proton = chi2;
					cdNDF_proton = NDF;
				}
			}
			HipoDataBank rectraj = (HipoDataBank) event.getBank("REC::Traj");
			for(int k = 0; k < rectraj.rows(); k++)
			{
				short pindex = rectraj.getShort("pindex", k);
				byte detector = rectraj.getByte("detector", k);
				byte layer = rectraj.getByte("layer", k);
				float x = rectraj.getFloat("x", k);
				float y = rectraj.getFloat("y", k);
				float z = rectraj.getFloat("z", k);
				if(pindex == e_pindex && detector == 6 && e_detector == 6 && layer ==  6)
				{
					dc_hit_e = new Vector3D (x, y, z);
					theta_deg_e = (float) (57.3*dc_hit_e.theta());
					phi_deg_e = (float) (57.3*dc_hit_e.phi());
					if(phi_deg_e <= -150) phi_deg_e = phi_deg_e+360;
					
					double[] b0 = thetaCorFac.getItem(e_sector, 0);
					double[] b1 = thetaCorFac.getItem(e_sector, 1);
					
					double a0 = (b0[4]*theta_deg_e*theta_deg_e*theta_deg_e*theta_deg_e)+(b0[3]*theta_deg_e*theta_deg_e*theta_deg_e)
									+(b0[2]*theta_deg_e*theta_deg_e)+(b0[1]*theta_deg_e)+b0[0];
					double a1 = (b1[4]*theta_deg_e*theta_deg_e*theta_deg_e*theta_deg_e)+(b1[3]*theta_deg_e*theta_deg_e*theta_deg_e)
									+(b1[2]*theta_deg_e*theta_deg_e)+(b1[1]*theta_deg_e)+b1[0];
					
					dc_hit_rot_e = new Vector3D (x*Math.cos(phi_rot[e_sector-1]/57.3)+y*Math.sin(phi_rot[e_sector-1]/57.3),
							y*Math.cos(phi_rot[e_sector-1]/57.3)-x*Math.sin(phi_rot[e_sector-1]/57.3), z);
					phi_rot_deg_e = (float) (57.3*dc_hit_rot_e.phi());
					pcf = (a1*phi_rot_deg_e)+a0;
				}
			}
			
			
			if(ecount == 1 && protoncount == 1
					&& v_e.z() > -12 && v_e.z() < 7 && Math.abs(v_e.z()-v_proton.z()) < 2.5 + (2.5/p_proton.p())
					&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30))
					&& theta_deg_e > 7 && theta_deg_e <= 30)
				//	&& 57.3*p_proton.theta() < 75)
			{
				W = (float) Math.sqrt((M*M)+(2*M*(E-p_e.p()))+(2*E*p_e.p()*((p_e.pz()/p_e.p())-1)));
				W_p_cor = (float) Math.sqrt((M*M)+(2*M*(E-(p_e.p()*pcf)))+(2*E*p_e.p()*pcf*((p_e.pz()/p_e.p())-1)));
				
				helastic_theta_vs_p_e.fill(p_e.p(), 57.3*p_e.theta());
				helastic_theta_vs_phi_e.fill(57.3*p_e.phi(), 57.3*p_e.theta());
				hW.fill(W);
				hpcorW.fill(W_p_cor);
				histGroups_W_sec.getItem(e_sector-1).fill(W);
				histGroups_p_cor_W_sec.getItem(e_sector-1).fill(W_p_cor);
				double[] theta_bnd  = new double[] {7, 8.5, 10, 12, 14, 16, 18, 21, 24, 33};
				for(int thbin = 0; thbin < 9; thbin++)
				{
					if(theta_deg_e > theta_bnd[thbin] && theta_deg_e <= theta_bnd[thbin+1])
					{
						double[] phi_rot_bnd = new double[] {-30, -9, 0, 9, 30};
						for(int phbin = 0; phbin < 4; phbin++)
						{
							if(phi_rot_deg_e > phi_rot_bnd[phbin] && phi_rot_deg_e <= phi_rot_bnd[phbin+1])
							{
								histGroups_W_theta_phi_bin.getItem(e_sector-1, thbin, phbin).fill(W);
								histGroups_p_cor_W_theta_phi_bin.getItem(e_sector-1, thbin, phbin).fill(W_p_cor);
								break;
							}
						}
						break;
					}
				}
			}
		}
	}
	
	static void processEventpiplusn(DataEvent event, int corEventCounter, IndexedList<double[]> thetaCorFac) {
		
		if(event.hasBank("REC::Particle") && event.hasBank("REC::Track") && event.hasBank("REC::Traj"))
		{
			HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
			int ecount = 0;
			int e_pindex = -1;
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0);
			LorentzVector cor_p_e = new LorentzVector(0, 0, 0, 0);
			Vector3D v_e = new Vector3D (0, 0, 0);
			float theta_deg_e = 0;
			float phi_deg_e = -200;
			float phi_rot_deg_e = -200;
			double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
			int pipluscount = 0;
			int piplus_pindex = -1;
			LorentzVector p_piplus = new LorentzVector(0, 0, 0, 0);
			Vector3D v_piplus = new Vector3D (0, 0, 0);
			Vector3D dc_hit_e = new Vector3D (0, 0, 0);
			Vector3D dc_hit_rot_e = new Vector3D (0, 0, 0);
			int neutroncount = 0;
			int protoncount = 0;
			double pcf = 0;
			for(int i = 0; i < rec.rows(); i++)
			{
				int pid = rec.getInt("pid", i);
				float px = rec.getFloat("px", i);
				float py = rec.getFloat("py", i);
				float pz = rec.getFloat("pz", i);
				float vx = rec.getFloat("vx", i);
				float vy = rec.getFloat("vy", i);
				float vz = rec.getFloat("vz", i);
				float p = (float) Math.sqrt((px*px)+(py*py)+(pz*pz));
				if(pid == 11)
				{
					ecount++;
					if(ecount == 1)
					{
						p_e = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.000511*0.000511)));
						v_e = new Vector3D (vx, vy, vz);
						e_pindex = i;
					}
				}
				if(pid == 211)
				{
					pipluscount++;
					if(pipluscount == 1)
					{
						p_piplus = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.139570*0.139570)));
						v_piplus = new Vector3D (vx, vy, vz);
						piplus_pindex = i;
					}
				}
				if(pid == 2112)
				{
					neutroncount++;
				}
				if(pid == 2212)
				{
					protoncount++;
				}
			}
			HipoDataBank rectrac = (HipoDataBank) event.getBank("REC::Track");
			byte e_detector = 0;
			byte e_sector = 0;
			for(int j = 0; j < rectrac.rows(); j++)
			{
				short pindex = rectrac.getShort("pindex", j);
				byte detector = rectrac.getByte("detector", j);
				byte sector = rectrac.getByte("sector", j);
				if(pindex == e_pindex)
				{
					e_detector = detector;
					e_sector = sector;
				}
			}
			HipoDataBank rectraj = (HipoDataBank) event.getBank("REC::Traj");
			for(int k = 0; k < rectraj.rows(); k++)
			{
				short pindex = rectraj.getShort("pindex", k);
				byte detector = rectraj.getByte("detector", k);
				byte layer = rectraj.getByte("layer", k);
				float x = rectraj.getFloat("x", k);
				float y = rectraj.getFloat("y", k);
				float z = rectraj.getFloat("z", k);
				if(pindex == e_pindex && detector == 6 && e_detector == 6 && layer ==  6)
				{
					dc_hit_e = new Vector3D (x, y, z);
					theta_deg_e = (float) (57.3*dc_hit_e.theta());
					phi_deg_e = (float) (57.3*dc_hit_e.phi());
					if(phi_deg_e <= -150) phi_deg_e = phi_deg_e+360;
					
					double[] b0 = thetaCorFac.getItem(e_sector, 0);
					double[] b1 = thetaCorFac.getItem(e_sector, 1);
					
					double a0 = (b0[4]*theta_deg_e*theta_deg_e*theta_deg_e*theta_deg_e)+(b0[3]*theta_deg_e*theta_deg_e*theta_deg_e)
									+(b0[2]*theta_deg_e*theta_deg_e)+(b0[1]*theta_deg_e)+b0[0];
					double a1 = (b1[4]*theta_deg_e*theta_deg_e*theta_deg_e*theta_deg_e)+(b1[3]*theta_deg_e*theta_deg_e*theta_deg_e)
									+(b1[2]*theta_deg_e*theta_deg_e)+(b1[1]*theta_deg_e)+b1[0];
					
					dc_hit_rot_e = new Vector3D (x*Math.cos(phi_rot[e_sector-1]/57.3)+y*Math.sin(phi_rot[e_sector-1]/57.3),
							y*Math.cos(phi_rot[e_sector-1]/57.3)-x*Math.sin(phi_rot[e_sector-1]/57.3), z);
					phi_rot_deg_e = (float) (57.3*dc_hit_rot_e.phi());
					pcf = (a1*phi_rot_deg_e)+a0;
					cor_p_e = new LorentzVector(pcf*p_e.px(), pcf*p_e.py(), pcf*p_e.pz(),
													Math.sqrt((pcf*p_e.p()*pcf*p_e.p())+(0.000511*0.000511)));
				}
			}
			double E = 6.535;
			LorentzVector p_B = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			LorentzVector p_T = new LorentzVector(0, 0, 0, 0.938272);
			LorentzVector p_X_epiplus = new LorentzVector(0, 0, 0, 0);
			p_X_epiplus.add(p_B);
			p_X_epiplus.add(p_T);
			p_X_epiplus.sub(p_e);
			p_X_epiplus.sub(p_piplus);
			LorentzVector p_cor_p_X_epiplus = new LorentzVector(0, 0, 0, 0);
			p_cor_p_X_epiplus.add(p_B);
			p_cor_p_X_epiplus.add(p_T);
			p_cor_p_X_epiplus.sub(cor_p_e);
			p_cor_p_X_epiplus.sub(p_piplus);
			
			if(ecount == 1 && pipluscount == 1// && neutroncount == 1
				//	&& 57.3*p_e.theta()<-10.7*p_e.p()+71.4// && rec.rows() == 3
					&& v_e.z() > -12 && v_e.z() < 7 && v_piplus.z() > -12 && v_piplus.z() < 7
					&& theta_deg_e > 7 && theta_deg_e <= 30)
			{
				hpiplusn_theta_vs_p_e.fill(p_e.p(), 57.3*p_e.theta());
				hpiplusn_theta_vs_phi_e.fill(57.3*p_e.phi(), 57.3*p_e.theta());
				hmissing_n.fill(p_X_epiplus.mass());
				hpcor_missing_n.fill(p_cor_p_X_epiplus.mass());
				histGroups_missing_n_sec.getItem(e_sector-1).fill(p_X_epiplus.mass());
				histGroups_p_cor_missing_n_sec.getItem(e_sector-1).fill(p_cor_p_X_epiplus.mass());
				double[] theta_bnd  = new double[] {7, 8.5, 10, 12, 14, 16, 18, 21, 24, 33};
				for(int thbin = 0; thbin < 9; thbin++)
				{
					if(theta_deg_e > theta_bnd[thbin] && theta_deg_e <= theta_bnd[thbin+1])
					{
						double[] phi_rot_bnd = new double[] {-30, -9, 0, 9, 30};
						for(int phbin = 0; phbin < 4; phbin++)
						{
							if(phi_rot_deg_e > phi_rot_bnd[phbin] && phi_rot_deg_e <= phi_rot_bnd[phbin+1])
							{
								histGroups_missing_n_theta_phi_bin.getItem(e_sector-1, thbin, phbin).fill(p_X_epiplus.mass());
								histGroups_p_cor_missing_n_theta_phi_bin.getItem(e_sector-1, thbin, phbin).fill(p_cor_p_X_epiplus.mass());
								break;
							}
						}
						break;
					}
				}
				if(57.3*p_e.theta()<-10*p_e.p()+60)
				{
					hpiplusn_theta_vs_p_e_psreg1.fill(p_e.p(), 57.3*p_e.theta());
					hpiplusn_theta_vs_phi_e_psreg1.fill(57.3*p_e.phi(), 57.3*p_e.theta());
					hmissing_n_psreg1.fill(p_X_epiplus.mass());
					hpcor_missing_n_psreg1.fill(p_cor_p_X_epiplus.mass());
					histGroups_missing_n_psreg1_sec.getItem(e_sector-1).fill(p_X_epiplus.mass());
					histGroups_p_cor_missing_n_psreg1_sec.getItem(e_sector-1).fill(p_cor_p_X_epiplus.mass());
					for(int thbin = 0; thbin < 9; thbin++)
					{
						if(theta_deg_e > theta_bnd[thbin] && theta_deg_e <= theta_bnd[thbin+1])
						{
							double[] phi_rot_bnd = new double[] {-30, -9, 0, 9, 30};
							for(int phbin = 0; phbin < 4; phbin++)
							{
								if(phi_rot_deg_e > phi_rot_bnd[phbin] && phi_rot_deg_e <= phi_rot_bnd[phbin+1])
								{
									histGroups_missing_n_theta_phi_bin_psreg1.getItem(e_sector-1, thbin, phbin).fill(p_X_epiplus.mass());
									histGroups_p_cor_missing_n_theta_phi_bin_psreg1.getItem(e_sector-1, thbin, phbin).fill(p_cor_p_X_epiplus.mass());
									break;
								}
							}
							break;
						}
					}
				}
				if(57.3*p_e.theta()>=-10*p_e.p()+60)
				{
					hpiplusn_theta_vs_p_e_psreg2.fill(p_e.p(), 57.3*p_e.theta());
					hpiplusn_theta_vs_phi_e_psreg2.fill(57.3*p_e.phi(), 57.3*p_e.theta());
					hmissing_n_psreg2.fill(p_X_epiplus.mass());
					hpcor_missing_n_psreg2.fill(p_cor_p_X_epiplus.mass());
					histGroups_missing_n_psreg2_sec.getItem(e_sector-1).fill(p_X_epiplus.mass());
					histGroups_p_cor_missing_n_psreg2_sec.getItem(e_sector-1).fill(p_cor_p_X_epiplus.mass());
					for(int thbin = 0; thbin < 9; thbin++)
					{
						if(theta_deg_e > theta_bnd[thbin] && theta_deg_e <= theta_bnd[thbin+1])
						{
							double[] phi_rot_bnd = new double[] {-30, -9, 0, 9, 30};
							for(int phbin = 0; phbin < 4; phbin++)
							{
								if(phi_rot_deg_e > phi_rot_bnd[phbin] && phi_rot_deg_e <= phi_rot_bnd[phbin+1])
								{
									histGroups_missing_n_theta_phi_bin_psreg2.getItem(e_sector-1, thbin, phbin).fill(p_X_epiplus.mass());
									histGroups_p_cor_missing_n_theta_phi_bin_psreg2.getItem(e_sector-1, thbin, phbin).fill(p_cor_p_X_epiplus.mass());
									break;
								}
							}
							break;
						}
					}
				}
			}
		}
	}

	static void processEventdvcs(DataEvent event, int corEventCounter, IndexedList<double[]> thetaCorFac) {
	
	if(event.hasBank("REC::Particle") && event.hasBank("REC::Track") && event.hasBank("REC::Traj"))
	{
		HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
		int ecount = 0;
		int e_pindex = -1;
		LorentzVector p_e = new LorentzVector(0, 0, 0, 0);
		LorentzVector cor_p_e = new LorentzVector(0, 0, 0, 0);
		Vector3D v_e = new Vector3D (0, 0, 0);
		float theta_deg_e = 0;
		float phi_deg_e = -200;
		float phi_rot_deg_e = -200;
		double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
		int protoncount = 0;
		int proton_pindex = -1;
		LorentzVector p_proton = new LorentzVector(0, 0, 0, 0);
		Vector3D v_proton = new Vector3D (0, 0, 0);
		float cdchi2_proton = 100;
		short cdNDF_proton = 100;
		int gammacount = 0;
		int gamma_pindex = -1;		
		LorentzVector p_gamma = new LorentzVector(0, 0, 0, 0);
		double cX = 0;
		double cY = 0;
		boolean fid_ecal_gamma = false;
		double[] cX_split  = new double[] {87, 82, 85, 77, 78, 82};
		double[] t_left  = new double[] {58.7356, 62.8204, 62.2296, 53.7756, 58.2888, 54.5822};
		double[] t_right  = new double[] {58.7477, 51.2589, 59.2357, 56.2415, 60.8219, 49.8914};
		double[] s_left  = new double[] {0.582053, 0.544976, 0.549788, 0.56899, 0.56414, 0.57343};
		double[] s_right  = new double[] {-0.591876, -0.562926, -0.562246, -0.563726, -0.568902, -0.550729};
		double[] r_left  = new double[] {64.9348, 64.7541, 67.832, 55.9324, 55.9225, 60.0997};
		double[] r_right  = new double[] {65.424, 54.6992, 63.6628, 57.8931, 56.5367, 56.4641};
		double[] q_left  = new double[] {0.745578, 0.606081, 0.729202, 0.627239, 0.503674, 0.717899};
		double[] q_right  = new double[] {-0.775022, -0.633863, -0.678901, -0.612458, -0.455319, -0.692481};
		Vector3D dc_hit_e = new Vector3D (0, 0, 0);
		Vector3D dc_hit_rot_e = new Vector3D (0, 0, 0);
		double pcf = 0;
		for(int i = 0; i < rec.rows(); i++)
		{
			int pid = rec.getInt("pid", i);
			float px = rec.getFloat("px", i);
			float py = rec.getFloat("py", i);
			float pz = rec.getFloat("pz", i);
			float vx = rec.getFloat("vx", i);
			float vy = rec.getFloat("vy", i);
			float vz = rec.getFloat("vz", i);
			float p = (float) Math.sqrt((px*px)+(py*py)+(pz*pz));
			if(pid == 11)
			{
				ecount++;
				if(ecount == 1)
				{
					p_e = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.000511*0.000511)));
					v_e = new Vector3D (vx, vy, vz);
					e_pindex = i;
				}
			}
			if(pid == 2212)
			{
				protoncount++;
				if(protoncount == 1)
				{
					p_proton = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.938272*0.938272)));
					v_proton = new Vector3D (vx, vy, vz);
					proton_pindex = i;
				}
			}
			if(pid == 22)
			{
				gammacount++;
				p_gamma = new LorentzVector(px, py, pz, p);
				gamma_pindex = i;
				
				//Start: gamma fiducial cut
				if(event.hasBank("REC::Calorimeter"))
				{
					HipoDataBank reccal = (HipoDataBank) event.getBank("REC::Calorimeter");
					for(int l = 0; l < reccal.rows(); l++)
					{
						short pindex = reccal.getShort("pindex", l);
						byte detector = reccal.getByte("detector", l);
						byte layer = reccal.getByte("layer", l);
						float x = reccal.getFloat("x", l);
						float y = reccal.getFloat("y", l);
						float lu = reccal.getFloat("lu", l);
						float lv = reccal.getFloat("lv", l);
						float lw = reccal.getFloat("lw", l);
						byte sector = reccal.getByte("sector", l);
						if(pindex == gamma_pindex && detector == 7 && layer == 1)
						{
							cX = x*Math.cos(rc*phi_rot[sector-1])+y*Math.sin(rc*phi_rot[sector-1]);
							cY = y*Math.cos(rc*phi_rot[sector-1])-x*Math.sin(rc*phi_rot[sector-1]);
							if(cX > cX_split[sector-1] && cY < s_left[sector-1]*(cX-t_left[sector-1])
								&& cY > s_right[sector-1]*(cX-t_right[sector-1]))
								fid_ecal_gamma = true;
							else if(cX < cX_split[sector-1] && cY < q_left[sector-1]*(cX-r_left[sector-1])
										&& cY > q_right[sector-1]*(cX-r_right[sector-1]))
								fid_ecal_gamma = true;
						}
						
						if(pindex == gamma_pindex && detector == 7 && sector == 1)
						{
							if(layer == 1)
							{
								if((lw > 74.5 && lw < 79.5) || (lw > 85.5 && lw < 90.5) || (lw > 213 && lw < 218)
										|| (lw > 224.5 && lw < 229.5) || (lw > 410)) fid_ecal_gamma = false;
							}
							else if(layer == 4)
							{
								if((lv > 70 && lv < 93)) fid_ecal_gamma = false;
							}
							else if(layer == 7)
							{
								if((lu > 410.5)) fid_ecal_gamma = false;
								if((lv < 24) || (lv > 36.5 && lv < 60) || (lv > 411)) fid_ecal_gamma = false;
								if((lw < 21.5)) fid_ecal_gamma = false;
							}
						}
						else if(pindex == gamma_pindex && detector == 7 && sector == 2)
						{
							if(layer == 1)
							{
								if((lv > 102 && lv < 113)) fid_ecal_gamma = false;
							}
							else if(layer == 4)
							{
								if((lu > 396)) fid_ecal_gamma = false;
								if((lw > 363)) fid_ecal_gamma = false;
							}
							else if(layer == 7)
							{
								if((lu < 12)) fid_ecal_gamma = false;
								if((lw < 10.5) || (lw > 376)) fid_ecal_gamma = false;
							}
						}
						else if(pindex == gamma_pindex && detector == 7 && sector == 3)
						{
							if(layer == 1)
							{
								if((lv > 354.5 && lv < 376.5)) fid_ecal_gamma = false;
							}
							else if(layer == 4)
							{
								if((lu < 23)) fid_ecal_gamma = false;
								if((lw < 10) || (lw > 363)) fid_ecal_gamma = false;
							}
							else if(layer == 7)
							{
								if((lw > 387)) fid_ecal_gamma = false;
							}
						}
						else if(pindex == gamma_pindex && detector == 7 && sector == 4)
						{
							if(layer == 1)
							{
								if((lv < 13) || (lv > 230 && lv < 240.5)) fid_ecal_gamma = false;
								if((lw > 410)) fid_ecal_gamma = false;
							}
							else if(layer == 4)
							{
								if((lu < 20.5)) fid_ecal_gamma = false;
							}
							else if(layer == 7)
							{
								if((lw < 32.5)) fid_ecal_gamma = false;
							}
						}
						else if(pindex == gamma_pindex && detector == 7 && sector == 5)
						{
							if(layer == 4)
							{
								if((lv < 23)) fid_ecal_gamma = false;
								if((lw < 10)) fid_ecal_gamma = false;
							}
							else if(layer == 7)
							{
								if((lu > 193.5 && lu < 217)) fid_ecal_gamma = false;
								if((lv < 24)) fid_ecal_gamma = false;
							}
						}
						else if(pindex == gamma_pindex && detector == 7 && sector == 6)
						{
							if(layer == 1)
							{
								if((lw > 174.5 && lw < 179.5) || (lw > 185.5 && lw < 190.5)) fid_ecal_gamma = false;
							}
							else if(layer == 4)
							{
								if((lv < 11.5)) fid_ecal_gamma = false;
								if((lu < 20.5)) fid_ecal_gamma = false;
							}
							else if(layer == 7)
							{
								if((lv < 12) || (lv > 423)) fid_ecal_gamma = false;
								if((lw < 32.5)) fid_ecal_gamma = false;
							}
						}		
					}
				}
				// End: gamma fiducial cut
			}	
		}
		HipoDataBank rectrac = (HipoDataBank) event.getBank("REC::Track");
		byte e_detector = 0;
		byte e_sector = 0;
		byte proton_detector = 0;
		for(int j = 0; j < rectrac.rows(); j++)
		{
			short pindex = rectrac.getShort("pindex", j);
			byte detector = rectrac.getByte("detector", j);
			byte sector = rectrac.getByte("sector", j);
			float chi2 = rectrac.getFloat("chi2", j);
			short NDF = rectrac.getShort("NDF", j);
			if(pindex == e_pindex)
			{
				e_detector = detector;
				e_sector = sector;
			}
			if(pindex == proton_pindex)
			{
				proton_detector = detector;
				cdchi2_proton = chi2;
				cdNDF_proton = NDF;
			}
		}
		HipoDataBank rectraj = (HipoDataBank) event.getBank("REC::Traj");
		for(int k = 0; k < rectraj.rows(); k++)
		{
			short pindex = rectraj.getShort("pindex", k);
			byte detector = rectraj.getByte("detector", k);
			byte layer = rectraj.getByte("layer", k);
			float x = rectraj.getFloat("x", k);
			float y = rectraj.getFloat("y", k);
			float z = rectraj.getFloat("z", k);
			if(pindex == e_pindex && detector == 6 && e_detector == 6 && layer ==  6)
			{
				dc_hit_e = new Vector3D (x, y, z);
				theta_deg_e = (float) (57.3*dc_hit_e.theta());
				phi_deg_e = (float) (57.3*dc_hit_e.phi());
				if(phi_deg_e <= -150) phi_deg_e = phi_deg_e+360;
				
				double[] b0 = thetaCorFac.getItem(e_sector, 0);
				double[] b1 = thetaCorFac.getItem(e_sector, 1);
				
				double a0 = (b0[4]*theta_deg_e*theta_deg_e*theta_deg_e*theta_deg_e)+(b0[3]*theta_deg_e*theta_deg_e*theta_deg_e)
								+(b0[2]*theta_deg_e*theta_deg_e)+(b0[1]*theta_deg_e)+b0[0];
				double a1 = (b1[4]*theta_deg_e*theta_deg_e*theta_deg_e*theta_deg_e)+(b1[3]*theta_deg_e*theta_deg_e*theta_deg_e)
								+(b1[2]*theta_deg_e*theta_deg_e)+(b1[1]*theta_deg_e)+b1[0];
				
				dc_hit_rot_e = new Vector3D (x*Math.cos(phi_rot[e_sector-1]/57.3)+y*Math.sin(phi_rot[e_sector-1]/57.3),
												y*Math.cos(phi_rot[e_sector-1]/57.3)-x*Math.sin(phi_rot[e_sector-1]/57.3), z);
				phi_rot_deg_e = (float) (57.3*dc_hit_rot_e.phi());
				pcf = (a1*phi_rot_deg_e)+a0;
				cor_p_e = new LorentzVector(pcf*p_e.px(), pcf*p_e.py(), pcf*p_e.pz(),
												Math.sqrt((pcf*p_e.p()*pcf*p_e.p())+(0.000511*0.000511)));
			}
		}
		
		
		double E = 6.535;
		LorentzVector p_B = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
		LorentzVector p_T = new LorentzVector(0, 0, 0, 0.938272);
		
		LorentzVector p_X_ep = new LorentzVector(0, 0, 0, 0);
		p_X_ep.add(p_B);
		p_X_ep.add(p_T);
		p_X_ep.sub(p_e);
		p_X_ep.sub(p_proton);
		LorentzVector p_X_eg = new LorentzVector(0, 0, 0, 0);
		p_X_eg.add(p_B);
		p_X_eg.add(p_T);
		p_X_eg.sub(p_e);
		p_X_eg.sub(p_gamma);
		LorentzVector p_X_epg = new LorentzVector(0, 0, 0, 0);
		p_X_epg.add(p_B);
		p_X_epg.add(p_T);
		p_X_epg.sub(p_e);
		p_X_epg.sub(p_proton);
		p_X_epg.sub(p_gamma);
		LorentzVector q = new LorentzVector(0, 0, 0, 0);
		q.add(p_B);
		q.sub(p_e);
		double Q2 = -1*q.mass2();
		LorentzVector w = new LorentzVector(0, 0, 0, 0);
		w.add(p_B);
		w.add(p_T);
		w.sub(p_e);
		double W = w.mass();
		LorentzVector p_cor_p_X_ep = new LorentzVector(0, 0, 0, 0);
		p_cor_p_X_ep.add(p_B);
		p_cor_p_X_ep.add(p_T);
		p_cor_p_X_ep.sub(cor_p_e);
		p_cor_p_X_ep.sub(p_proton);
		LorentzVector p_cor_p_X_eg = new LorentzVector(0, 0, 0, 0);
		p_cor_p_X_eg.add(p_B);
		p_cor_p_X_eg.add(p_T);
		p_cor_p_X_eg.sub(cor_p_e);
		p_cor_p_X_eg.sub(p_gamma);
		LorentzVector p_cor_p_X_epg = new LorentzVector(0, 0, 0, 0);
		p_cor_p_X_epg.add(p_B);
		p_cor_p_X_epg.add(p_T);
		p_cor_p_X_epg.sub(cor_p_e);
		p_cor_p_X_epg.sub(p_proton);
		p_cor_p_X_epg.sub(p_gamma);
		LorentzVector p_cor_q = new LorentzVector(0, 0, 0, 0);
		p_cor_q.add(p_B);
		p_cor_q.sub(cor_p_e);
		double p_cor_Q2 = -1*p_cor_q.mass2();
		LorentzVector p_cor_w = new LorentzVector(0, 0, 0, 0);
		p_cor_w.add(p_B);
		p_cor_w.add(p_T);
		p_cor_w.sub(cor_p_e);
		double p_cor_W = p_cor_w.mass();
		
		if(ecount == 1 && protoncount == 1 && gammacount == 1
				&& v_e.z() > -12 && v_e.z() < 7 && Math.abs(v_e.z()-v_proton.z()) < 2.5 + (2.5/p_proton.p())
				&& p_gamma.vect().theta(p_e.vect()) > 5
				&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30))
				&& 57.3*p_proton.theta() < 75
				&& fid_ecal_gamma == true
				&& theta_deg_e > 7 && theta_deg_e <= 30)
		{
			if(Q2 > 1 && W > 2)
			{
				 	hdvcs_theta_vs_p_e.fill(p_e.p(), 57.3*p_e.theta());
				 	hdvcs_theta_vs_phi_e.fill( 57.3*p_e.phi(), 57.3*p_e.theta());
					htheta_cone_gamma.fill(p_gamma.vect().theta(p_X_ep.vect()));
					hX_eg_m.fill(p_X_eg.mass());
					hX_epg_pt.fill(Math.sqrt(p_X_epg.px()*p_X_epg.px()+p_X_epg.py()*p_X_epg.py()));
					hX_epg_E.fill(p_X_epg.e());
					histGroups_theta_cone_gamma_sec.getItem(e_sector-1).fill(p_gamma.vect().theta(p_X_ep.vect()));				
					histGroups_X_eg_m_sec.getItem(e_sector-1).fill(p_X_eg.mass());				
					histGroups_X_epg_pt_sec .getItem(e_sector-1).fill(Math.sqrt(p_X_epg.px()*p_X_epg.px()+p_X_epg.py()*p_X_epg.py()));				
					histGroups_X_epg_E_sec.getItem(e_sector-1).fill(p_X_epg.e());
			}
			if(p_cor_Q2 > 1 && p_cor_W > 2)
			{
					hpcor_theta_cone_gamma.fill(p_gamma.vect().theta(p_cor_p_X_ep.vect()));
					hpcor_X_eg_m.fill(p_cor_p_X_eg.mass());
					hpcor_X_epg_pt.fill(Math.sqrt(p_cor_p_X_epg.px()*p_cor_p_X_epg.px()+p_cor_p_X_epg.py()*p_cor_p_X_epg.py()));
					hpcor_X_epg_E.fill(p_cor_p_X_epg.e());
					histGroups_p_cor_theta_cone_gamma_sec.getItem(e_sector-1).fill(p_gamma.vect().theta(p_cor_p_X_ep.vect()));
					histGroups_p_cor_X_eg_m_sec.getItem(e_sector-1).fill(p_cor_p_X_eg.mass());
					histGroups_p_cor_X_epg_pt_sec.getItem(e_sector-1).fill(Math.sqrt(p_cor_p_X_epg.px()*p_cor_p_X_epg.px()
																						+p_cor_p_X_epg.py()*p_cor_p_X_epg.py()));
					histGroups_p_cor_X_epg_E_sec.getItem(e_sector-1).fill(p_cor_p_X_epg.e());
			}
			double[] theta_bnd  = new double[] {7, 8.5, 10, 12, 14, 16, 18, 21, 24, 33};
			for(int thbin = 0; thbin < 9; thbin++)
			{
				if(theta_deg_e > theta_bnd[thbin] && theta_deg_e <= theta_bnd[thbin+1])
				{
					double[] phi_rot_bnd = new double[] {-30, -9, 0, 9, 30};
					for(int phbin = 0; phbin < 4; phbin++)
					{
						if(phi_rot_deg_e > phi_rot_bnd[phbin] && phi_rot_deg_e <= phi_rot_bnd[phbin+1])
						{
							if(Q2 > 1 && W > 2)
							{
								histGroups_theta_cone_gamma_theta_phi_bin.getItem(e_sector-1, thbin, phbin)
																			.fill(p_gamma.vect().theta(p_X_ep.vect()));	
								histGroups_X_eg_m_theta_phi_bin.getItem(e_sector-1, thbin, phbin).fill(p_X_eg.mass());	
								histGroups_X_epg_pt_theta_phi_bin.getItem(e_sector-1, thbin, phbin)
																	.fill(Math.sqrt(p_X_epg.px()*p_X_epg.px()+p_X_epg.py()*p_X_epg.py()));		
								histGroups_X_epg_E_theta_phi_bin.getItem(e_sector-1, thbin, phbin).fill(p_X_epg.e());
							}
							if(p_cor_Q2 > 1 && p_cor_W > 2)
							{
								histGroups_p_cor_theta_cone_gamma_theta_phi_bin.getItem(e_sector-1, thbin, phbin)
																				.fill(p_gamma.vect().theta(p_cor_p_X_ep.vect()));
								histGroups_p_cor_X_eg_m_theta_phi_bin.getItem(e_sector-1, thbin, phbin).fill(p_cor_p_X_eg.mass());
								histGroups_p_cor_X_epg_pt_theta_phi_bin.getItem(e_sector-1, thbin, phbin)
																		.fill(Math.sqrt(p_cor_p_X_epg.px()*p_cor_p_X_epg.px()
																						+p_cor_p_X_epg.py()*p_cor_p_X_epg.py()));
								histGroups_p_cor_X_epg_E_theta_phi_bin.getItem(e_sector-1, thbin, phbin).fill(p_cor_p_X_epg.e());
							}
							break;
						}
					}
					break;
				}
			}
		}
	}
}

	public static void main(String[] args) {
		
		histos_W_sec();
		histos_W_vs_phi_bin();
		histos_missing_n_sec();
		histos_missing_n_psreg1_sec();
		histos_missing_n_psreg2_sec();
		histos_missing_n_vs_phi_bin();
		histos_missing_n_vs_phi_bin_psreg1();
		histos_missing_n_vs_phi_bin_psreg2();
		histos_dvcs_ev_sec();
		histos_dvcs_ev_vs_phi_bin();

		IndexedList<double[]> thetaCorFac = new IndexedList<>(2);
		double[] thetaCorFac_sec1_b0 = new double[] {0.9938077481295297, 8.258708800790782E-4, -2.2658532616476987E-5,
														9.677115140791733E-7, -4.5642462769384737E-8};
		double[] thetaCorFac_sec1_b1 = new double[] {0.0024178549687052706, -6.653055638986705E-4, 7.219407837197144E-5,
														-3.180462703382774E-6, 4.9352098479104826E-8};
		double[] thetaCorFac_sec2_b0 = new double[] {1.0035661097100197, -0.003285642462705683, 5.292556662659571E-4,
														-3.0231148934779212E-5, 5.609421915775323E-7};
		double[] thetaCorFac_sec2_b1 = new double[] {5.346964410100895E-4, -6.714148268661577E-5, 5.201475751662365E-6,
														-1.6155087929089417E-7, 1.0186220422620223E-9};
		double[] thetaCorFac_sec3_b0 = new double[] {1.0087429705449467, -0.004408627402638188, 6.482030554858202E-4,
														-3.394142844711124E-5, 5.779315008382153E-7};
		double[] thetaCorFac_sec3_b1 = new double[] {-4.397136504767579E-4, 1.4884197765913248E-4, -1.8433837666625867E-5,
														8.244944205916959E-7, -1.228942642526823E-8};
		double[] thetaCorFac_sec4_b0 = new double[] {1.006956918896026, -0.003711867838701133, 5.518086095685475E-4,
														-2.8746618476848303E-5, 4.828563256249679E-7};
		double[] thetaCorFac_sec4_b1 = new double[] {8.309626628189099E-4, -2.701812621232045E-4, 2.875524544466852E-5,
														-1.2835027495862537E-6, 2.029581677685505E-8};
		double[] thetaCorFac_sec5_b0 = new double[] {0.9570390799595814, 0.009328675535642626, -6.714343815803992E-4,
														1.987029852607302E-5, -2.1755105045221692E-7};
		double[] thetaCorFac_sec5_b1 = new double[] {-5.873762378695643E-4, 1.9121282318863288E-4, -1.9101091097978772E-5,
														7.82569493538631E-7, -1.1798923261435443E-8};
		double[] thetaCorFac_sec6_b0 = new double[] {0.9636208140401962, 0.00948121762243233, -7.335233524854341E-4,
														2.4169978390757133E-5, -3.071125516271284E-7};
		double[] thetaCorFac_sec6_b1 = new double[] {0.0017450577912736217, -5.245913063775676E-4, 5.3940149514332456E-5,
														-2.3107384729861627E-6, 3.575938540124231E-8};
		thetaCorFac.add(thetaCorFac_sec1_b0, 1, 0);
		thetaCorFac.add(thetaCorFac_sec1_b1, 1, 1);
		thetaCorFac.add(thetaCorFac_sec2_b0, 2, 0);
		thetaCorFac.add(thetaCorFac_sec2_b1, 2, 1);
		thetaCorFac.add(thetaCorFac_sec3_b0, 3, 0);
		thetaCorFac.add(thetaCorFac_sec3_b1, 3, 1);
		thetaCorFac.add(thetaCorFac_sec4_b0, 4, 0);
		thetaCorFac.add(thetaCorFac_sec4_b1, 4, 1);
		thetaCorFac.add(thetaCorFac_sec5_b0, 5, 0);
		thetaCorFac.add(thetaCorFac_sec5_b1, 5, 1);
		thetaCorFac.add(thetaCorFac_sec6_b0, 6, 0);
		thetaCorFac.add(thetaCorFac_sec6_b1, 6, 1);
		
		double[] theta_bnd  = new double[] {7, 8.5, 10, 12, 14, 16, 18, 21, 24, 33};
		double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
		double[] phirotba = new double[] {-13.27, -4.05, 4.05, 13.27};
		/*
	// Elastic	
		
		readerelastic.open("C:/Users/joshtanj/Documents/download/merged_skim_ep_bank_6535MeV_skim_elastic_59xx.hipo");
		
		int elasticEventCounter = 0;
		while(readerelastic.hasEvent())// && elasticEventCounter < 10000000)
		{
			elasticEventCounter++;
			processEventelastic(readerelastic.getNextEvent(), elasticEventCounter, thetaCorFac);
			if(elasticEventCounter%1000000 == 0) System.out.println("Elastic Event: " + elasticEventCounter);
		}

		readerelastic.close();
		
		F1D fW;
		
		JFrame frameelastickin = new JFrame("Elastic All Sector Kinematics");
		frameelastickin.setSize(1000, 1000);
		EmbeddedCanvas canelastickin = new EmbeddedCanvas();
		frameelastickin.add(canelastickin);
		frameelastickin.setLocationRelativeTo(null);
		frameelastickin.setVisible(true);
		canelastickin.divide(2, 2);
		canelastickin.cd(0);
		canelastickin.setFont("Arial");
		helastic_theta_vs_p_e.setTitle("e^- #theta vs. p");
		helastic_theta_vs_p_e.setTitleX("p [GeV]");
		helastic_theta_vs_p_e.setTitleY("#theta [#degree]");
		canelastickin.getPad(0).setTitleFontSize(32);
		canelastickin.getPad(0).setAxisTitleFontSize(32);
		canelastickin.getPad(0).setAxisLabelFontSize(24);
		canelastickin.draw(helastic_theta_vs_p_e, "same");
		canelastickin.cd(1);
		canelastickin.setFont("Arial");
		helastic_theta_vs_phi_e.setTitle("e^- #theta vs. #phi");
		helastic_theta_vs_phi_e.setTitleX("#phi [#degree]");
		helastic_theta_vs_phi_e.setTitleY("#theta [#degree]");
		canelastickin.getPad(1).setTitleFontSize(32);
		canelastickin.getPad(1).setAxisTitleFontSize(32);
		canelastickin.getPad(1).setAxisLabelFontSize(24);
		canelastickin.draw(helastic_theta_vs_phi_e, "same");
		canelastickin.cd(2);
		canelastickin.setFont("Arial");
		hW.setTitle("All Sector W");
		hW.setTitleX("W [GeV]");
		hW.setTitleY("Counts");
		hW.setOptStat(10);
		canelastickin.getPad(2).setTitleFontSize(32);
		canelastickin.getPad(2).setAxisTitleFontSize(32);
		canelastickin.getPad(2).setAxisLabelFontSize(24);
		canelastickin.getPad(2).setStatBoxFontSize(18);
		canelastickin.draw(hW, "same");
		fW = new F1D("fW", "[amp]*gaus(x,[mean],[sigma])",
							0.7+(hW.getMaximumBin()*(1.2-0.7)/100)-0.035,
							0.7+(hW.getMaximumBin()*(1.2-0.7)/100)+0.025);
		fW.setParameter(0, hW.getMax()*1.2);
		fW.setParameter(1, 0.7+(hW.getMaximumBin()*(1.2-0.7)/100));
		fW.setParameter(2, 0.01);
		DataFitter.fit(fW, hW, "Q");
		fW.setLineColor(2);
		fW.setLineWidth(3);
		fW.setOptStat(11110);
		canelastickin.draw(fW, "same");
		canelastickin.cd(3);
		canelastickin.setFont("Arial");
		hpcorW.setTitle("All Sector p_e Corrected W");
		hpcorW.setTitleX("W [GeV]");
		hpcorW.setTitleY("Counts");
		hpcorW.setOptStat(10);
		canelastickin.getPad(3).setTitleFontSize(32);
		canelastickin.getPad(3).setAxisTitleFontSize(32);
		canelastickin.getPad(3).setAxisLabelFontSize(24);
		canelastickin.getPad(3).setStatBoxFontSize(18);
		canelastickin.draw(hpcorW, "same");
		fW = new F1D("fW", "[amp]*gaus(x,[mean],[sigma])",
						0.7+(hpcorW.getMaximumBin()*(1.2-0.7)/100)-0.035,
						0.7+(hpcorW.getMaximumBin()*(1.2-0.7)/100)+0.025);
		fW.setParameter(0, hpcorW.getMax()*1.2);
		fW.setParameter(1, 0.7+(hpcorW.getMaximumBin()*(1.2-0.7)/100));
		fW.setParameter(2, 0.01);
		DataFitter.fit(fW, hpcorW, "Q");
		fW.setLineColor(2);
		fW.setLineWidth(3);
		fW.setOptStat(11110);
		canelastickin.draw(fW, "same");
		
		JFrame frameWsec = new JFrame("W");
		frameWsec.setSize(1500, 1000);
		EmbeddedCanvas canWsec = new EmbeddedCanvas();
		frameWsec.add(canWsec);
		frameWsec.setLocationRelativeTo(null);
		frameWsec.setVisible(true);
		canWsec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canWsec.cd(seci);
			canWsec.setFont("Arial");
			histGroups_W_sec.getItem(seci).setTitle("Sector " + (seci+1) + " W");
			histGroups_W_sec.getItem(seci).setTitleX("W [GeV]");
			histGroups_W_sec.getItem(seci).setTitleY("Counts");
			histGroups_W_sec.getItem(seci).setOptStat(10);
			canWsec.getPad(seci).setTitleFontSize(32);
			canWsec.getPad(seci).setAxisTitleFontSize(32);
			canWsec.getPad(seci).setAxisLabelFontSize(24);
			canWsec.getPad(seci).setStatBoxFontSize(18);
			canWsec.draw(histGroups_W_sec.getItem(seci), "same");
			fW = new F1D("fW", "[amp]*gaus(x,[mean],[sigma])",
								0.7+(histGroups_W_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)-0.035,
								0.7+(histGroups_W_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)+0.025);
			fW.setParameter(0, histGroups_W_sec.getItem(seci).getMax()*1.2);
			fW.setParameter(1, 0.7+(histGroups_W_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100));
			fW.setParameter(2, 0.01);
			DataFitter.fit(fW, histGroups_W_sec.getItem(seci), "Q");
			fW.setLineColor(2);
			fW.setLineWidth(3);
			fW.setOptStat(11110);
			canWsec.draw(fW, "same");
		}
		
		JFrame framepcorWsec = new JFrame("p_e Corrected W");
		framepcorWsec.setSize(1500, 1000);
		EmbeddedCanvas canpcorWsec = new EmbeddedCanvas();
		framepcorWsec.add(canpcorWsec);
		framepcorWsec.setLocationRelativeTo(null);
		framepcorWsec.setVisible(true);
		canpcorWsec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canpcorWsec.cd(seci);
			canpcorWsec.setFont("Arial");
			histGroups_p_cor_W_sec.getItem(seci).setTitle("Sector " + (seci+1) + " W");
			histGroups_p_cor_W_sec.getItem(seci).setTitleX("W [GeV]");
			histGroups_p_cor_W_sec.getItem(seci).setTitleY("Counts");
			histGroups_p_cor_W_sec.getItem(seci).setOptStat(10);
			canpcorWsec.getPad(seci).setTitleFontSize(32);
			canpcorWsec.getPad(seci).setAxisTitleFontSize(32);
			canpcorWsec.getPad(seci).setAxisLabelFontSize(24);
			canpcorWsec.getPad(seci).setStatBoxFontSize(18);
			canpcorWsec.draw(histGroups_p_cor_W_sec.getItem(seci), "same");
			fW = new F1D("fW", "[amp]*gaus(x,[mean],[sigma])",
								0.7+(histGroups_p_cor_W_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)-0.035,
								0.7+(histGroups_p_cor_W_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)+0.025);
			fW.setParameter(0, histGroups_p_cor_W_sec.getItem(seci).getMax()*1.2);
			fW.setParameter(1, 0.7+(histGroups_p_cor_W_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100));
			fW.setParameter(2, 0.01);
			DataFitter.fit(fW, histGroups_p_cor_W_sec.getItem(seci), "Q");
			fW.setLineColor(2);
			fW.setLineWidth(3);
			fW.setOptStat(11110);
			canpcorWsec.draw(fW, "same");
		}
		
		IndexedList<Double> meanW = new IndexedList<>(3);
		IndexedList<Double> sigmaW = new IndexedList<>(3);
		
		IndexedList<Double> meanpcorW = new IndexedList<>(3);
		IndexedList<Double> sigmapcorW = new IndexedList<>(3);
		
		for(int seci = 0; seci < 6; seci++)
		{
			for(int theta_bini = 0; theta_bini < 9; theta_bini++)
			{
				for(int phi_bini = 0; phi_bini < 4; phi_bini++)
				{
					fW = new F1D("fW", "[amp]*gaus(x,[mean],[sigma])", 
										0.7+(histGroups_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()
												*(1.2-0.7)/100)-0.035,
										0.7+(histGroups_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()
												*(1.2-0.7)/100)+0.025);
					fW.setParameter(0, histGroups_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax());
					fW.setParameter(1, 0.7+(histGroups_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()
												*(1.2-0.7)/100));
					fW.setParameter(2, 0.01);
					DataFitter.fit(fW, histGroups_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini), "Q");
					meanW.add(fW.getParameter(1), seci, theta_bini, phi_bini);
					sigmaW.add(fW.getParameter(2), seci, theta_bini, phi_bini);
					fW = new F1D("fW", "[amp]*gaus(x,[mean],[sigma])", 
										0.7+(histGroups_p_cor_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()
												*(1.2-0.7)/100)-0.035,
										0.7+(histGroups_p_cor_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()
												*(1.2-0.7)/100)+0.025);
					fW.setParameter(0, histGroups_p_cor_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax());
					fW.setParameter(1, 0.7+(histGroups_p_cor_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()
										*(1.2-0.7)/100));
					fW.setParameter(2, 0.01);
					DataFitter.fit(fW, histGroups_p_cor_W_theta_phi_bin.getItem(seci, theta_bini, phi_bini), "Q");
					meanpcorW.add(fW.getParameter(1), seci, theta_bini, phi_bini);
					sigmapcorW.add(fW.getParameter(2), seci, theta_bini, phi_bini);
				}
			}
		}
		
		IndexedList<GraphErrors> histerrGroups_W_vs_phi = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> histerrGroups_p_cor_W_vs_phi = new IndexedList<GraphErrors>(2);
		
		for(int seci = 0; seci < 6; seci++) {
			for(int theta_bini = 0; theta_bini < 9; theta_bini++) {
				GraphErrors hW_vs_phi = new GraphErrors();
				histerrGroups_W_vs_phi.add(hW_vs_phi, seci, theta_bini);
				GraphErrors hp_cor_W_vs_phi = new GraphErrors();
				histerrGroups_p_cor_W_vs_phi.add(hp_cor_W_vs_phi, seci, theta_bini);
				for(int phi_bini = 0; phi_bini < 4; phi_bini++)
				{
					if(meanW.hasItem(seci, theta_bini,  phi_bini) && sigmaW.hasItem(seci, theta_bini, phi_bini)
							&& sigmaW.getItem(seci, theta_bini, phi_bini) < 0.2)
						histerrGroups_W_vs_phi.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]), meanW.getItem(seci, theta_bini, phi_bini),
															0, sigmaW.getItem(seci, theta_bini, phi_bini));
					if(meanpcorW.hasItem(seci, theta_bini,  phi_bini) && sigmapcorW.hasItem(seci, theta_bini, phi_bini)
							&& sigmapcorW.getItem(seci, theta_bini, phi_bini) < 0.2)
						histerrGroups_p_cor_W_vs_phi.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]), meanpcorW.getItem(seci, theta_bini, phi_bini),
															0, sigmapcorW.getItem(seci, theta_bini, phi_bini));
				}
			}
		}
		
		for(int theta_bini = 0; theta_bini < 9; theta_bini++){
			JFrame frameWphibin = new JFrame("W vs. #phi Binned");
			frameWphibin.setSize(1500, 1000);
			EmbeddedCanvas canWphibin = new EmbeddedCanvas();
			frameWphibin.add(canWphibin);
			frameWphibin.setLocationRelativeTo(null);
			frameWphibin.setVisible(true);
			canWphibin.divide(3, 2);
			for(int seci = 0; seci < 6; seci++){
				canWphibin.cd(seci);
				canWphibin.setFont("Arial");
				histerrGroups_W_vs_phi.getItem(seci, theta_bini).setTitle("e^- Hit: " + theta_bnd[theta_bini] + "#degree < #theta <= "
																		+ theta_bnd[theta_bini+1] + " #degree, W vs. #phi");
				histerrGroups_W_vs_phi.getItem(seci, theta_bini).setTitleX("#phi [#degree]");
				histerrGroups_W_vs_phi.getItem(seci, theta_bini).setTitleY("Sector " + (seci+1) + " W: #mu +/- #sigma [GeV]");
				histerrGroups_W_vs_phi.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_W_vs_phi.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_W_vs_phi.getItem(seci, theta_bini).setLineThickness(1);
				canWphibin.getPad(seci).setTitleFontSize(32);
				canWphibin.getPad(seci).setAxisTitleFontSize(32);
				canWphibin.getPad(seci).setAxisLabelFontSize(24);
				canWphibin.getPad(seci).setAxisRange((phi_rot[seci]-30), (phi_rot[seci]+30), 0.7, 1.2);
				histerrGroups_p_cor_W_vs_phi.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_p_cor_W_vs_phi.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_p_cor_W_vs_phi.getItem(seci, theta_bini).setMarkerColor(2);
				histerrGroups_p_cor_W_vs_phi.getItem(seci, theta_bini).setLineColor(2);
				histerrGroups_p_cor_W_vs_phi.getItem(seci, theta_bini).setLineThickness(1);
				canWphibin.draw(histerrGroups_W_vs_phi.getItem(seci, theta_bini), "same");
				canWphibin.draw(histerrGroups_p_cor_W_vs_phi.getItem(seci, theta_bini), "same");
			} 
		}
		*/
		/*
	// pi+n	
		
		readerpiplusn.open("C:/Users/joshtanj/Documents/download/merged_skim_epi+_bank_6535MeV_skim_custom_5885.hipo");
		
		int piplusnEventCounter = 0;
		while(readerpiplusn.hasEvent() && piplusnEventCounter < 20000000)
		{
			piplusnEventCounter++;
			processEventpiplusn(readerpiplusn.getNextEvent(), piplusnEventCounter, thetaCorFac);
			if(piplusnEventCounter%5000000 == 0) System.out.println("ep->s#pi^+n Event: " + piplusnEventCounter);
		}

		readerpiplusn.close();
		
		F1D fnm;
		
		JFrame framepiplusnkin = new JFrame("ep->e#pi^+n All Sector Kinematics");
		framepiplusnkin.setSize(1000, 1000);
		EmbeddedCanvas canpiplusnkin = new EmbeddedCanvas();
		framepiplusnkin.add(canpiplusnkin);
		framepiplusnkin.setLocationRelativeTo(null);
		framepiplusnkin.setVisible(true);
		canpiplusnkin.divide(2, 2);
		canpiplusnkin.cd(0);
		canpiplusnkin.setFont("Arial");
		hpiplusn_theta_vs_p_e.setTitle("e^- #theta vs. p");
		hpiplusn_theta_vs_p_e.setTitleX("p [GeV]");
		hpiplusn_theta_vs_p_e.setTitleY("#theta [#degree]");
		canpiplusnkin.getPad(0).setTitleFontSize(32);
		canpiplusnkin.getPad(0).setAxisTitleFontSize(32);
		canpiplusnkin.getPad(0).setAxisLabelFontSize(24);
		canpiplusnkin.draw(hpiplusn_theta_vs_p_e, "same");
		canpiplusnkin.cd(1);
		canpiplusnkin.setFont("Arial");
		hpiplusn_theta_vs_phi_e.setTitle("e^- #theta vs. #phi");
		hpiplusn_theta_vs_phi_e.setTitleX("#phi [#degree]");
		hpiplusn_theta_vs_phi_e.setTitleY("#theta [#degree]");
		canpiplusnkin.getPad(1).setTitleFontSize(32);
		canpiplusnkin.getPad(1).setAxisTitleFontSize(32);
		canpiplusnkin.getPad(1).setAxisLabelFontSize(24);
		canpiplusnkin.draw(hpiplusn_theta_vs_phi_e, "same");
		canpiplusnkin.cd(2);
		canpiplusnkin.setFont("Arial");
		hmissing_n.setTitle("All Sector n_miss Mass");
		hmissing_n.setTitleX("M [GeV]");
		hmissing_n.setTitleY("Counts");
		hmissing_n.setOptStat(10);
		canpiplusnkin.getPad(2).setTitleFontSize(32);
		canpiplusnkin.getPad(2).setAxisTitleFontSize(32);
		canpiplusnkin.getPad(2).setAxisLabelFontSize(24);
		canpiplusnkin.getPad(2).setStatBoxFontSize(18);
		canpiplusnkin.draw(hmissing_n, "same");
		fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
								0.7+(hmissing_n.getMaximumBin()*(1.2-0.7)/100)-0.045,
								0.7+(hmissing_n.getMaximumBin()*(1.2-0.7)/100)+0.04);
		fnm.setParameter(0, hmissing_n.getMax()*1.2);
		fnm.setParameter(1, 0.7+(hmissing_n.getMaximumBin()*(1.2-0.7)/100));
		fnm.setParameter(2, 0.03);
		DataFitter.fit(fnm, hmissing_n, "Q");
		fnm.setLineColor(2);
		fnm.setLineWidth(3);
		fnm.setOptStat(11110);
		canpiplusnkin.draw(fnm, "same");
		canpiplusnkin.cd(3);
		canpiplusnkin.setFont("Arial");
		hpcor_missing_n.setTitle("All Sector p_e Corrected n_miss Mass");
		hpcor_missing_n.setTitleX("M [GeV]");
		hpcor_missing_n.setTitleY("Counts");
		hpcor_missing_n.setOptStat(10);
		canpiplusnkin.getPad(3).setTitleFontSize(32);
		canpiplusnkin.getPad(3).setAxisTitleFontSize(32);
		canpiplusnkin.getPad(3).setAxisLabelFontSize(24);
		canpiplusnkin.getPad(3).setStatBoxFontSize(18);
		canpiplusnkin.draw(hpcor_missing_n, "same");
		fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
						0.7+(hpcor_missing_n.getMaximumBin()*(1.2-0.7)/100)-0.045,
						0.7+(hpcor_missing_n.getMaximumBin()*(1.2-0.7)/100)+0.04);
		fnm.setParameter(0, hpcor_missing_n.getMax()*1.2);
		fnm.setParameter(1, 0.7+(hpcor_missing_n.getMaximumBin()*(1.2-0.7)/100));
		fnm.setParameter(2, 0.03);
		DataFitter.fit(fnm, hpcor_missing_n, "Q");
		fnm.setLineColor(2);
		fnm.setLineWidth(3);
		fnm.setOptStat(11110);
		canpiplusnkin.draw(fnm, "same");
		/*
		JFrame framemissingnsec = new JFrame("n_miss Mass");
		framemissingnsec.setSize(1500, 1000);
		EmbeddedCanvas canmissingnsec = new EmbeddedCanvas();
		framemissingnsec.add(canmissingnsec);
		framemissingnsec.setLocationRelativeTo(null);
		framemissingnsec.setVisible(true);
		canmissingnsec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canmissingnsec.cd(seci);
			canmissingnsec.setFont("Arial");
			histGroups_missing_n_sec.getItem(seci).setTitle("Sector " + (seci+1) + " n_miss Mass");
			histGroups_missing_n_sec.getItem(seci).setTitleX("M [GeV]");
			histGroups_missing_n_sec.getItem(seci).setTitleY("Counts");
			histGroups_missing_n_sec.getItem(seci).setOptStat(10);
			canmissingnsec.getPad(seci).setTitleFontSize(32);
			canmissingnsec.getPad(seci).setAxisTitleFontSize(32);
			canmissingnsec.getPad(seci).setAxisLabelFontSize(24);
			canmissingnsec.getPad(seci).setStatBoxFontSize(18);
			canmissingnsec.draw(histGroups_missing_n_sec.getItem(seci), "same");
			fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
									0.7+(histGroups_missing_n_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)-0.045,
									0.7+(histGroups_missing_n_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)+0.04);
			fnm.setParameter(0, histGroups_missing_n_sec.getItem(seci).getMax()*1.2);
			fnm.setParameter(1, 0.7+(histGroups_missing_n_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100));
			fnm.setParameter(2, 0.03);
			DataFitter.fit(fnm, histGroups_missing_n_sec.getItem(seci), "Q");
			fnm.setLineColor(2);
			fnm.setLineWidth(3);
			fnm.setOptStat(11110);
			canmissingnsec.draw(fnm, "same");
		}
		
		JFrame framepcormissingnsec = new JFrame("p_e Corrected n_miss Mass");
		framepcormissingnsec.setSize(1500, 1000);
		EmbeddedCanvas canpcormissingnsec = new EmbeddedCanvas();
		framepcormissingnsec.add(canpcormissingnsec);
		framepcormissingnsec.setLocationRelativeTo(null);
		framepcormissingnsec.setVisible(true);
		canpcormissingnsec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canpcormissingnsec.cd(seci);
			canpcormissingnsec.setFont("Arial");
			histGroups_p_cor_missing_n_sec.getItem(seci).setTitle("Sector " + (seci+1) + " n_miss Mass");
			histGroups_p_cor_missing_n_sec.getItem(seci).setTitleX("M [GeV]");
			histGroups_p_cor_missing_n_sec.getItem(seci).setTitleY("Counts");
			histGroups_p_cor_missing_n_sec.getItem(seci).setOptStat(10);
			canpcormissingnsec.getPad(seci).setTitleFontSize(32);
			canpcormissingnsec.getPad(seci).setAxisTitleFontSize(32);
			canpcormissingnsec.getPad(seci).setAxisLabelFontSize(24);
			canpcormissingnsec.getPad(seci).setStatBoxFontSize(18);
			canpcormissingnsec.draw(histGroups_p_cor_missing_n_sec.getItem(seci), "same");
			fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
									0.7+(histGroups_p_cor_missing_n_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)-0.045,
									0.7+(histGroups_p_cor_missing_n_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)+0.04);
			fnm.setParameter(0, histGroups_p_cor_missing_n_sec.getItem(seci).getMax()*1.2);
			fnm.setParameter(1, 0.7+(histGroups_p_cor_missing_n_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100));
			fnm.setParameter(2, 0.03);
			DataFitter.fit(fnm, histGroups_p_cor_missing_n_sec.getItem(seci), "Q");
			fnm.setLineColor(2);
			fnm.setLineWidth(3);
			fnm.setOptStat(11110);
			canpcormissingnsec.draw(fnm, "same");
		}
		
		IndexedList<Double> meanmissingn = new IndexedList<>(3);
		IndexedList<Double> sigmamissingn = new IndexedList<>(3);
		
		IndexedList<Double> meanpcormissingn = new IndexedList<>(3);
		IndexedList<Double> sigmapcormissingn = new IndexedList<>(3);
		
		for(int seci = 0; seci < 6; seci++)
		{
			for(int theta_bini = 0; theta_bini < 9; theta_bini++)
			{
				for(int phi_bini = 0; phi_bini < 4; phi_bini++)
				{
					fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
											0.7+(histGroups_missing_n_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()
													*(1.2-0.7)/100)-0.045,
											0.7+(histGroups_missing_n_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()
													*(1.2-0.7)/100)+0.04);
					fnm.setParameter(0, histGroups_missing_n_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax()*1.2);
					fnm.setParameter(1, 0.7+(histGroups_missing_n_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()
										*(1.2-0.7)/100));
					fnm.setParameter(2, 0.03);
					DataFitter.fit(fnm, histGroups_missing_n_theta_phi_bin.getItem(seci, theta_bini, phi_bini), "Q");
					meanmissingn.add(fnm.getParameter(1), seci, theta_bini, phi_bini);
					sigmamissingn.add(fnm.getParameter(2), seci, theta_bini, phi_bini);
					fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
											0.7+(histGroups_p_cor_missing_n_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()
													*(1.2-0.7)/100)-0.045,
											0.7+(histGroups_p_cor_missing_n_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()
													*(1.2-0.7)/100)+0.04);
					fnm.setParameter(0, histGroups_p_cor_missing_n_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax()*1.2);
					fnm.setParameter(1, 0.7+(histGroups_p_cor_missing_n_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()
										*(1.2-0.7)/100));
					fnm.setParameter(2, 0.03);
					DataFitter.fit(fnm, histGroups_p_cor_missing_n_theta_phi_bin.getItem(seci, theta_bini, phi_bini), "Q");
					meanpcormissingn.add(fnm.getParameter(1), seci, theta_bini, phi_bini);
					sigmapcormissingn.add(fnm.getParameter(2), seci, theta_bini, phi_bini);
				}
			}
		}
		
		IndexedList<GraphErrors> histerrGroups_missing_n_vs_phi = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> histerrGroups_p_cor_missing_n_vs_phi = new IndexedList<GraphErrors>(2);
		
		for(int seci = 0; seci < 6; seci++) {
			for(int theta_bini = 0; theta_bini < 9; theta_bini++) {
				GraphErrors hmissing_n_vs_phi = new GraphErrors();
				histerrGroups_missing_n_vs_phi.add(hmissing_n_vs_phi, seci, theta_bini);
				GraphErrors hp_cor_missing_n_vs_phi = new GraphErrors();
				histerrGroups_p_cor_missing_n_vs_phi.add(hp_cor_missing_n_vs_phi, seci, theta_bini);
				for(int phi_bini = 0; phi_bini < 4; phi_bini++)
				{
					if(meanmissingn.hasItem(seci, theta_bini,  phi_bini) && sigmamissingn.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_missing_n_vs_phi.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]),
															meanmissingn.getItem(seci, theta_bini,  phi_bini),
															0, sigmamissingn.getItem(seci, theta_bini, phi_bini));
					if(meanpcormissingn.hasItem(seci, theta_bini,  phi_bini) && sigmapcormissingn.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_p_cor_missing_n_vs_phi.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]),
															meanpcormissingn.getItem(seci, theta_bini,  phi_bini),
															0, sigmapcormissingn.getItem(seci, theta_bini, phi_bini));
				}
			}
		}
		
		for(int theta_bini = 0; theta_bini < 9; theta_bini++){
			JFrame framemissingnphibin = new JFrame("n_miss Mass vs. #phi Binned");
			framemissingnphibin.setSize(1500, 1000);
			EmbeddedCanvas canmissingnphibin = new EmbeddedCanvas();
			framemissingnphibin.add(canmissingnphibin);
			framemissingnphibin.setLocationRelativeTo(null);
			framemissingnphibin.setVisible(true);
			canmissingnphibin.divide(3, 2);
			for(int seci = 0; seci < 6; seci++){
				canmissingnphibin.cd(seci);
				canmissingnphibin.setFont("Arial");
				histerrGroups_missing_n_vs_phi.getItem(seci, theta_bini).setTitle("e^- Hit: " + theta_bnd[theta_bini] + "#degree < #theta <= "
																		+ theta_bnd[theta_bini+1] + " #degree, n_miss Mass vs. #phi");
				histerrGroups_missing_n_vs_phi.getItem(seci, theta_bini).setTitleX("#phi [#degree]");
				histerrGroups_missing_n_vs_phi.getItem(seci, theta_bini).setTitleY("Sector " + (seci+1) + " n_miss : #mu +/- #sigma [GeV]");
				histerrGroups_missing_n_vs_phi.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_missing_n_vs_phi.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_missing_n_vs_phi.getItem(seci, theta_bini).setLineThickness(1);
				canmissingnphibin.getPad(seci).setTitleFontSize(32);
				canmissingnphibin.getPad(seci).setAxisTitleFontSize(32);
				canmissingnphibin.getPad(seci).setAxisLabelFontSize(24);
				canmissingnphibin.getPad(seci).setAxisRange((phi_rot[seci]-30), (phi_rot[seci]+30), 0.7, 1.2);
				histerrGroups_p_cor_missing_n_vs_phi.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_p_cor_missing_n_vs_phi.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_p_cor_missing_n_vs_phi.getItem(seci, theta_bini).setMarkerColor(2);
				histerrGroups_p_cor_missing_n_vs_phi.getItem(seci, theta_bini).setLineColor(2);
				histerrGroups_p_cor_missing_n_vs_phi.getItem(seci, theta_bini).setLineThickness(1);			
				canmissingnphibin.draw(histerrGroups_missing_n_vs_phi.getItem(seci, theta_bini), "same");
				canmissingnphibin.draw(histerrGroups_p_cor_missing_n_vs_phi.getItem(seci, theta_bini), "same");
			}
		}
		
		JFrame framepiplusnkinpsreg1 = new JFrame("#theta<-10p+60: ep->e#pi^+n All Sector Kinematics");
		framepiplusnkinpsreg1.setSize(1000, 1000);
		EmbeddedCanvas canpiplusnkinpsreg1 = new EmbeddedCanvas();
		framepiplusnkinpsreg1.add(canpiplusnkinpsreg1);
		framepiplusnkinpsreg1.setLocationRelativeTo(null);
		framepiplusnkinpsreg1.setVisible(true);
		canpiplusnkinpsreg1.divide(2, 2);
		canpiplusnkinpsreg1.cd(0);
		canpiplusnkinpsreg1.setFont("Arial");
		hpiplusn_theta_vs_p_e_psreg1.setTitle("e^- #theta vs. p");
		hpiplusn_theta_vs_p_e_psreg1.setTitleX("p [GeV]");
		hpiplusn_theta_vs_p_e_psreg1.setTitleY("#theta [#degree]");
		canpiplusnkinpsreg1.getPad(0).setTitleFontSize(32);
		canpiplusnkinpsreg1.getPad(0).setAxisTitleFontSize(32);
		canpiplusnkinpsreg1.getPad(0).setAxisLabelFontSize(24);
		canpiplusnkinpsreg1.draw(hpiplusn_theta_vs_p_e_psreg1, "same");
		canpiplusnkinpsreg1.cd(1);
		canpiplusnkinpsreg1.setFont("Arial");
		hpiplusn_theta_vs_phi_e.setTitle("e^- #theta vs. #phi");
		hpiplusn_theta_vs_phi_e.setTitleX("#phi [#degree]");
		hpiplusn_theta_vs_phi_e.setTitleY("#theta [#degree]");
		canpiplusnkinpsreg1.getPad(1).setTitleFontSize(32);
		canpiplusnkinpsreg1.getPad(1).setAxisTitleFontSize(32);
		canpiplusnkinpsreg1.getPad(1).setAxisLabelFontSize(24);
		canpiplusnkinpsreg1.draw(hpiplusn_theta_vs_phi_e, "same");
		canpiplusnkinpsreg1.cd(2);
		canpiplusnkinpsreg1.setFont("Arial");
		hmissing_n_psreg1.setTitle("All Sector n_miss Mass");
		hmissing_n_psreg1.setTitleX("M [GeV]");
		hmissing_n_psreg1.setTitleY("Counts");
		hmissing_n_psreg1.setOptStat(10);
		canpiplusnkinpsreg1.getPad(2).setTitleFontSize(32);
		canpiplusnkinpsreg1.getPad(2).setAxisTitleFontSize(32);
		canpiplusnkinpsreg1.getPad(2).setAxisLabelFontSize(24);
		canpiplusnkinpsreg1.getPad(2).setStatBoxFontSize(18);
		canpiplusnkinpsreg1.draw(hmissing_n_psreg1, "same");
		fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
								0.7+(hmissing_n_psreg1.getMaximumBin()*(1.2-0.7)/100)-0.045,
								0.7+(hmissing_n_psreg1.getMaximumBin()*(1.2-0.7)/100)+0.04);
		fnm.setParameter(0, hmissing_n_psreg1.getMax()*1.2);
		fnm.setParameter(1, 0.7+(hmissing_n_psreg1.getMaximumBin()*(1.2-0.7)/100));
		fnm.setParameter(2, 0.03);
		DataFitter.fit(fnm, hmissing_n_psreg1, "Q");
		fnm.setLineColor(2);
		fnm.setLineWidth(3);
		fnm.setOptStat(11110);
		canpiplusnkinpsreg1.draw(fnm, "same");
		canpiplusnkinpsreg1.cd(3);
		canpiplusnkinpsreg1.setFont("Arial");
		hpcor_missing_n_psreg1.setTitle("All Sector p_e Corrected n_miss Mass");
		hpcor_missing_n_psreg1.setTitleX("M [GeV]");
		hpcor_missing_n_psreg1.setTitleY("Counts");
		hpcor_missing_n_psreg1.setOptStat(10);
		canpiplusnkinpsreg1.getPad(3).setTitleFontSize(32);
		canpiplusnkinpsreg1.getPad(3).setAxisTitleFontSize(32);
		canpiplusnkinpsreg1.getPad(3).setAxisLabelFontSize(24);
		canpiplusnkinpsreg1.getPad(3).setStatBoxFontSize(18);
		canpiplusnkinpsreg1.draw(hpcor_missing_n_psreg1, "same");
		fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
						0.7+(hpcor_missing_n_psreg1.getMaximumBin()*(1.2-0.7)/100)-0.045,
						0.7+(hpcor_missing_n_psreg1.getMaximumBin()*(1.2-0.7)/100)+0.04);
		fnm.setParameter(0, hpcor_missing_n_psreg1.getMax()*1.2);
		fnm.setParameter(1, 0.7+(hpcor_missing_n_psreg1.getMaximumBin()*(1.2-0.7)/100));
		fnm.setParameter(2, 0.03);
		DataFitter.fit(fnm, hpcor_missing_n_psreg1, "Q");
		fnm.setLineColor(2);
		fnm.setLineWidth(3);
		fnm.setOptStat(11110);
		canpiplusnkinpsreg1.draw(fnm, "same");
		
		JFrame framemissingnpsreg1sec = new JFrame("#theta<-10p+60: n_miss Mass");
		framemissingnpsreg1sec.setSize(1500, 1000);
		EmbeddedCanvas canmissingnpsreg1sec = new EmbeddedCanvas();
		framemissingnpsreg1sec.add(canmissingnpsreg1sec);
		framemissingnpsreg1sec.setLocationRelativeTo(null);
		framemissingnpsreg1sec.setVisible(true);
		canmissingnpsreg1sec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canmissingnpsreg1sec.cd(seci);
			canmissingnpsreg1sec.setFont("Arial");
			histGroups_missing_n_psreg1_sec.getItem(seci).setTitle("Sector " + (seci+1) + " n_miss Mass");
			histGroups_missing_n_psreg1_sec.getItem(seci).setTitleX("M [GeV]");
			histGroups_missing_n_psreg1_sec.getItem(seci).setTitleY("Counts");
			histGroups_missing_n_psreg1_sec.getItem(seci).setOptStat(10);
			canmissingnpsreg1sec.getPad(seci).setTitleFontSize(32);
			canmissingnpsreg1sec.getPad(seci).setAxisTitleFontSize(32);
			canmissingnpsreg1sec.getPad(seci).setAxisLabelFontSize(24);
			canmissingnpsreg1sec.getPad(seci).setStatBoxFontSize(18);
			canmissingnpsreg1sec.draw(histGroups_missing_n_psreg1_sec.getItem(seci), "same");
			fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
									0.7+(histGroups_missing_n_psreg1_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)-0.045,
									0.7+(histGroups_missing_n_psreg1_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)+0.04);
			fnm.setParameter(0, histGroups_missing_n_psreg1_sec.getItem(seci).getMax()*1.2);
			fnm.setParameter(1, 0.7+(histGroups_missing_n_psreg1_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100));
			fnm.setParameter(2, 0.03);
			DataFitter.fit(fnm, histGroups_missing_n_psreg1_sec.getItem(seci), "Q");
			fnm.setLineColor(2);
			fnm.setLineWidth(3);
			fnm.setOptStat(11110);
			canmissingnpsreg1sec.draw(fnm, "same");
		}
		
		JFrame framepcormissingnpsreg1sec = new JFrame("#theta<-10p+60: p_e Corrected n_miss Mass");
		framepcormissingnpsreg1sec.setSize(1500, 1000);
		EmbeddedCanvas canpcormissingnpsreg1sec = new EmbeddedCanvas();
		framepcormissingnpsreg1sec.add(canpcormissingnpsreg1sec);
		framepcormissingnpsreg1sec.setLocationRelativeTo(null);
		framepcormissingnpsreg1sec.setVisible(true);
		canpcormissingnpsreg1sec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canpcormissingnpsreg1sec.cd(seci);
			canpcormissingnpsreg1sec.setFont("Arial");
			histGroups_p_cor_missing_n_psreg1_sec.getItem(seci).setTitle("Sector " + (seci+1) + " n_miss Mass");
			histGroups_p_cor_missing_n_psreg1_sec.getItem(seci).setTitleX("M [GeV]");
			histGroups_p_cor_missing_n_psreg1_sec.getItem(seci).setTitleY("Counts");
			histGroups_p_cor_missing_n_psreg1_sec.getItem(seci).setOptStat(10);
			canpcormissingnpsreg1sec.getPad(seci).setTitleFontSize(32);
			canpcormissingnpsreg1sec.getPad(seci).setAxisTitleFontSize(32);
			canpcormissingnpsreg1sec.getPad(seci).setAxisLabelFontSize(24);
			canpcormissingnpsreg1sec.getPad(seci).setStatBoxFontSize(18);
			canpcormissingnpsreg1sec.draw(histGroups_p_cor_missing_n_psreg1_sec.getItem(seci), "same");
			fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
									0.7+(histGroups_p_cor_missing_n_psreg1_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)-0.045,
									0.7+(histGroups_p_cor_missing_n_psreg1_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)+0.04);
			fnm.setParameter(0, histGroups_p_cor_missing_n_psreg1_sec.getItem(seci).getMax()*1.2);
			fnm.setParameter(1, 0.7+(histGroups_p_cor_missing_n_psreg1_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100));
			fnm.setParameter(2, 0.03);
			DataFitter.fit(fnm, histGroups_p_cor_missing_n_psreg1_sec.getItem(seci), "Q");
			fnm.setLineColor(2);
			fnm.setLineWidth(3);
			fnm.setOptStat(11110);
			canpcormissingnpsreg1sec.draw(fnm, "same");
		}
		
		IndexedList<Double> meanmissingnpsreg1 = new IndexedList<>(3);
		IndexedList<Double> sigmamissingnpsreg1 = new IndexedList<>(3);
		
		IndexedList<Double> meanpcormissingnpsreg1 = new IndexedList<>(3);
		IndexedList<Double> sigmapcormissingnpsreg1 = new IndexedList<>(3);
		
		for(int seci = 0; seci < 6; seci++)
		{
			for(int theta_bini = 0; theta_bini < 9; theta_bini++)
			{
				for(int phi_bini = 0; phi_bini < 4; phi_bini++)
				{
					fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
											0.7+(histGroups_missing_n_theta_phi_bin_psreg1.getItem(seci, theta_bini, phi_bini).getMaximumBin()
													*(1.2-0.7)/100)-0.045,
											0.7+(histGroups_missing_n_theta_phi_bin_psreg1.getItem(seci, theta_bini, phi_bini).getMaximumBin()
													*(1.2-0.7)/100)+0.04);
					fnm.setParameter(0, histGroups_missing_n_theta_phi_bin_psreg1.getItem(seci, theta_bini, phi_bini).getMax()*1.2);
					fnm.setParameter(1, 0.7+(histGroups_missing_n_theta_phi_bin_psreg1.getItem(seci, theta_bini, phi_bini).getMaximumBin()
										*(1.2-0.7)/100));
					fnm.setParameter(2, 0.03);
					DataFitter.fit(fnm, histGroups_missing_n_theta_phi_bin_psreg1.getItem(seci, theta_bini, phi_bini), "Q");
					meanmissingnpsreg1.add(fnm.getParameter(1), seci, theta_bini, phi_bini);
					sigmamissingnpsreg1.add(fnm.getParameter(2), seci, theta_bini, phi_bini);
					fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
											0.7+(histGroups_p_cor_missing_n_theta_phi_bin_psreg1.getItem(seci, theta_bini, phi_bini).getMaximumBin()
													*(1.2-0.7)/100)-0.045,
											0.7+(histGroups_p_cor_missing_n_theta_phi_bin_psreg1.getItem(seci, theta_bini, phi_bini).getMaximumBin()
													*(1.2-0.7)/100)+0.04);
					fnm.setParameter(0, histGroups_p_cor_missing_n_theta_phi_bin_psreg1.getItem(seci, theta_bini, phi_bini).getMax()*1.2);
					fnm.setParameter(1, 0.7+(histGroups_p_cor_missing_n_theta_phi_bin_psreg1.getItem(seci, theta_bini, phi_bini).getMaximumBin()
										*(1.2-0.7)/100));
					fnm.setParameter(2, 0.03);
					DataFitter.fit(fnm, histGroups_p_cor_missing_n_theta_phi_bin_psreg1.getItem(seci, theta_bini, phi_bini), "Q");
					meanpcormissingnpsreg1.add(fnm.getParameter(1), seci, theta_bini, phi_bini);
					sigmapcormissingnpsreg1.add(fnm.getParameter(2), seci, theta_bini, phi_bini);
				}
			}
		}
		
		IndexedList<GraphErrors> histerrGroups_missing_n_vs_phi_psreg1 = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> histerrGroups_p_cor_missing_n_vs_phi_psreg1 = new IndexedList<GraphErrors>(2);
		
		for(int seci = 0; seci < 6; seci++) {
			for(int theta_bini = 0; theta_bini < 9; theta_bini++) {
				GraphErrors hmissing_n_vs_phi_psreg1 = new GraphErrors();
				histerrGroups_missing_n_vs_phi_psreg1.add(hmissing_n_vs_phi_psreg1, seci, theta_bini);
				GraphErrors hp_cor_missing_n_vs_phi_psreg1 = new GraphErrors();
				histerrGroups_p_cor_missing_n_vs_phi_psreg1.add(hp_cor_missing_n_vs_phi_psreg1, seci, theta_bini);
				for(int phi_bini = 0; phi_bini < 4; phi_bini++)
				{
					if(meanmissingnpsreg1.hasItem(seci, theta_bini,  phi_bini) && sigmamissingnpsreg1.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_missing_n_vs_phi_psreg1.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]),
															meanmissingnpsreg1.getItem(seci, theta_bini,  phi_bini),
															0, sigmamissingnpsreg1.getItem(seci, theta_bini, phi_bini));
					if(meanpcormissingnpsreg1.hasItem(seci, theta_bini,  phi_bini) && sigmapcormissingnpsreg1.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_p_cor_missing_n_vs_phi_psreg1.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]),
															meanpcormissingnpsreg1.getItem(seci, theta_bini,  phi_bini),
															0, sigmapcormissingnpsreg1.getItem(seci, theta_bini, phi_bini));
				}
			}
		}
		
		for(int theta_bini = 0; theta_bini < 9; theta_bini++){
			JFrame framemissingnpsreg1phibin = new JFrame("#theta<-10p+60: n_miss Mass vs. #phi Binned");
			framemissingnpsreg1phibin.setSize(1500, 1000);
			EmbeddedCanvas canmissingnpsreg1phibin = new EmbeddedCanvas();
			framemissingnpsreg1phibin.add(canmissingnpsreg1phibin);
			framemissingnpsreg1phibin.setLocationRelativeTo(null);
			framemissingnpsreg1phibin.setVisible(true);
			canmissingnpsreg1phibin.divide(3, 2);
			for(int seci = 0; seci < 6; seci++){
				canmissingnpsreg1phibin.cd(seci);
				canmissingnpsreg1phibin.setFont("Arial");
				histerrGroups_missing_n_vs_phi_psreg1.getItem(seci, theta_bini).setTitle("e^- Hit: " + theta_bnd[theta_bini] + "#degree < #theta <= "
																		+ theta_bnd[theta_bini+1] + " #degree, n_miss Mass vs. #phi");
				histerrGroups_missing_n_vs_phi_psreg1.getItem(seci, theta_bini).setTitleX("#phi [#degree]");
				histerrGroups_missing_n_vs_phi_psreg1.getItem(seci, theta_bini).setTitleY("Sector " + (seci+1) + " n_miss : #mu +/- #sigma [GeV]");
				histerrGroups_missing_n_vs_phi_psreg1.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_missing_n_vs_phi_psreg1.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_missing_n_vs_phi_psreg1.getItem(seci, theta_bini).setLineThickness(1);
				canmissingnpsreg1phibin.getPad(seci).setTitleFontSize(32);
				canmissingnpsreg1phibin.getPad(seci).setAxisTitleFontSize(32);
				canmissingnpsreg1phibin.getPad(seci).setAxisLabelFontSize(24);
				canmissingnpsreg1phibin.getPad(seci).setAxisRange((phi_rot[seci]-30), (phi_rot[seci]+30), 0.7, 1.2);
				histerrGroups_p_cor_missing_n_vs_phi_psreg1.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_p_cor_missing_n_vs_phi_psreg1.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_p_cor_missing_n_vs_phi_psreg1.getItem(seci, theta_bini).setMarkerColor(2);
				histerrGroups_p_cor_missing_n_vs_phi_psreg1.getItem(seci, theta_bini).setLineColor(2);
				histerrGroups_p_cor_missing_n_vs_phi_psreg1.getItem(seci, theta_bini).setLineThickness(1);			
				canmissingnpsreg1phibin.draw(histerrGroups_missing_n_vs_phi_psreg1.getItem(seci, theta_bini), "same");
				canmissingnpsreg1phibin.draw(histerrGroups_p_cor_missing_n_vs_phi_psreg1.getItem(seci, theta_bini), "same");
			}
		}
		
		JFrame framepiplusnkinpsreg2 = new JFrame("#theta>=-10p+60: ep->e#pi^+n All Sector Kinematics");
		framepiplusnkinpsreg2.setSize(1000, 1000);
		EmbeddedCanvas canpiplusnkinpsreg2 = new EmbeddedCanvas();
		framepiplusnkinpsreg2.add(canpiplusnkinpsreg2);
		framepiplusnkinpsreg2.setLocationRelativeTo(null);
		framepiplusnkinpsreg2.setVisible(true);
		canpiplusnkinpsreg2.divide(2, 2);
		canpiplusnkinpsreg2.cd(0);
		canpiplusnkinpsreg2.setFont("Arial");
		hpiplusn_theta_vs_p_e_psreg2.setTitle("e^- #theta vs. p");
		hpiplusn_theta_vs_p_e_psreg2.setTitleX("p [GeV]");
		hpiplusn_theta_vs_p_e_psreg2.setTitleY("#theta [#degree]");
		canpiplusnkinpsreg2.getPad(0).setTitleFontSize(32);
		canpiplusnkinpsreg2.getPad(0).setAxisTitleFontSize(32);
		canpiplusnkinpsreg2.getPad(0).setAxisLabelFontSize(24);
		canpiplusnkinpsreg2.draw(hpiplusn_theta_vs_p_e_psreg2, "same");
		canpiplusnkinpsreg2.cd(1);
		canpiplusnkinpsreg2.setFont("Arial");
		hpiplusn_theta_vs_phi_e.setTitle("e^- #theta vs. #phi");
		hpiplusn_theta_vs_phi_e.setTitleX("#phi [#degree]");
		hpiplusn_theta_vs_phi_e.setTitleY("#theta [#degree]");
		canpiplusnkinpsreg2.getPad(1).setTitleFontSize(32);
		canpiplusnkinpsreg2.getPad(1).setAxisTitleFontSize(32);
		canpiplusnkinpsreg2.getPad(1).setAxisLabelFontSize(24);
		canpiplusnkinpsreg2.draw(hpiplusn_theta_vs_phi_e, "same");
		canpiplusnkinpsreg2.cd(2);
		canpiplusnkinpsreg2.setFont("Arial");
		hmissing_n_psreg2.setTitle("All Sector n_miss Mass");
		hmissing_n_psreg2.setTitleX("M [GeV]");
		hmissing_n_psreg2.setTitleY("Counts");
		hmissing_n_psreg2.setOptStat(10);
		canpiplusnkinpsreg2.getPad(2).setTitleFontSize(32);
		canpiplusnkinpsreg2.getPad(2).setAxisTitleFontSize(32);
		canpiplusnkinpsreg2.getPad(2).setAxisLabelFontSize(24);
		canpiplusnkinpsreg2.getPad(2).setStatBoxFontSize(18);
		canpiplusnkinpsreg2.draw(hmissing_n_psreg2, "same");
		fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
								0.7+(hmissing_n_psreg2.getMaximumBin()*(1.2-0.7)/100)-0.045,
								0.7+(hmissing_n_psreg2.getMaximumBin()*(1.2-0.7)/100)+0.04);
		fnm.setParameter(0, hmissing_n_psreg2.getMax()*1.2);
		fnm.setParameter(1, 0.7+(hmissing_n_psreg2.getMaximumBin()*(1.2-0.7)/100));
		fnm.setParameter(2, 0.03);
		DataFitter.fit(fnm, hmissing_n_psreg2, "Q");
		fnm.setLineColor(2);
		fnm.setLineWidth(3);
		fnm.setOptStat(11110);
		canpiplusnkinpsreg2.draw(fnm, "same");
		canpiplusnkinpsreg2.cd(3);
		canpiplusnkinpsreg2.setFont("Arial");
		hpcor_missing_n_psreg2.setTitle("All Sector p_e Corrected n_miss Mass");
		hpcor_missing_n_psreg2.setTitleX("M [GeV]");
		hpcor_missing_n_psreg2.setTitleY("Counts");
		hpcor_missing_n_psreg2.setOptStat(10);
		canpiplusnkinpsreg2.getPad(3).setTitleFontSize(32);
		canpiplusnkinpsreg2.getPad(3).setAxisTitleFontSize(32);
		canpiplusnkinpsreg2.getPad(3).setAxisLabelFontSize(24);
		canpiplusnkinpsreg2.getPad(3).setStatBoxFontSize(18);
		canpiplusnkinpsreg2.draw(hpcor_missing_n_psreg2, "same");
		fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
						0.7+(hpcor_missing_n_psreg2.getMaximumBin()*(1.2-0.7)/100)-0.045,
						0.7+(hpcor_missing_n_psreg2.getMaximumBin()*(1.2-0.7)/100)+0.04);
		fnm.setParameter(0, hpcor_missing_n_psreg2.getMax()*1.2);
		fnm.setParameter(1, 0.7+(hpcor_missing_n_psreg2.getMaximumBin()*(1.2-0.7)/100));
		fnm.setParameter(2, 0.03);
		DataFitter.fit(fnm, hpcor_missing_n_psreg2, "Q");
		fnm.setLineColor(2);
		fnm.setLineWidth(3);
		fnm.setOptStat(11110);
		canpiplusnkinpsreg2.draw(fnm, "same");
		
		JFrame framemissingnpsreg2sec = new JFrame("#theta>=-10p+60: n_miss Mass");
		framemissingnpsreg2sec.setSize(1500, 1000);
		EmbeddedCanvas canmissingnpsreg2sec = new EmbeddedCanvas();
		framemissingnpsreg2sec.add(canmissingnpsreg2sec);
		framemissingnpsreg2sec.setLocationRelativeTo(null);
		framemissingnpsreg2sec.setVisible(true);
		canmissingnpsreg2sec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canmissingnpsreg2sec.cd(seci);
			canmissingnpsreg2sec.setFont("Arial");
			histGroups_missing_n_psreg2_sec.getItem(seci).setTitle("Sector " + (seci+1) + " n_miss Mass");
			histGroups_missing_n_psreg2_sec.getItem(seci).setTitleX("M [GeV]");
			histGroups_missing_n_psreg2_sec.getItem(seci).setTitleY("Counts");
			histGroups_missing_n_psreg2_sec.getItem(seci).setOptStat(10);
			canmissingnpsreg2sec.getPad(seci).setTitleFontSize(32);
			canmissingnpsreg2sec.getPad(seci).setAxisTitleFontSize(32);
			canmissingnpsreg2sec.getPad(seci).setAxisLabelFontSize(24);
			canmissingnpsreg2sec.getPad(seci).setStatBoxFontSize(18);
			canmissingnpsreg2sec.draw(histGroups_missing_n_psreg2_sec.getItem(seci), "same");
			fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
									0.7+(histGroups_missing_n_psreg2_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)-0.045,
									0.7+(histGroups_missing_n_psreg2_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)+0.04);
			fnm.setParameter(0, histGroups_missing_n_psreg2_sec.getItem(seci).getMax()*1.2);
			fnm.setParameter(1, 0.7+(histGroups_missing_n_psreg2_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100));
			fnm.setParameter(2, 0.03);
			DataFitter.fit(fnm, histGroups_missing_n_psreg2_sec.getItem(seci), "Q");
			fnm.setLineColor(2);
			fnm.setLineWidth(3);
			fnm.setOptStat(11110);
			canmissingnpsreg2sec.draw(fnm, "same");
		}
		
		JFrame framepcormissingnpsreg2sec = new JFrame("#theta>=-10p+60: p_e Corrected n_miss Mass");
		framepcormissingnpsreg2sec.setSize(1500, 1000);
		EmbeddedCanvas canpcormissingnpsreg2sec = new EmbeddedCanvas();
		framepcormissingnpsreg2sec.add(canpcormissingnpsreg2sec);
		framepcormissingnpsreg2sec.setLocationRelativeTo(null);
		framepcormissingnpsreg2sec.setVisible(true);
		canpcormissingnpsreg2sec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canpcormissingnpsreg2sec.cd(seci);
			canpcormissingnpsreg2sec.setFont("Arial");
			histGroups_p_cor_missing_n_psreg2_sec.getItem(seci).setTitle("Sector " + (seci+1) + " n_miss Mass");
			histGroups_p_cor_missing_n_psreg2_sec.getItem(seci).setTitleX("M [GeV]");
			histGroups_p_cor_missing_n_psreg2_sec.getItem(seci).setTitleY("Counts");
			histGroups_p_cor_missing_n_psreg2_sec.getItem(seci).setOptStat(10);
			canpcormissingnpsreg2sec.getPad(seci).setTitleFontSize(32);
			canpcormissingnpsreg2sec.getPad(seci).setAxisTitleFontSize(32);
			canpcormissingnpsreg2sec.getPad(seci).setAxisLabelFontSize(24);
			canpcormissingnpsreg2sec.getPad(seci).setStatBoxFontSize(18);
			canpcormissingnpsreg2sec.draw(histGroups_p_cor_missing_n_psreg2_sec.getItem(seci), "same");
			fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
									0.7+(histGroups_p_cor_missing_n_psreg2_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)-0.045,
									0.7+(histGroups_p_cor_missing_n_psreg2_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100)+0.04);
			fnm.setParameter(0, histGroups_p_cor_missing_n_psreg2_sec.getItem(seci).getMax()*1.2);
			fnm.setParameter(1, 0.7+(histGroups_p_cor_missing_n_psreg2_sec.getItem(seci).getMaximumBin()*(1.2-0.7)/100));
			fnm.setParameter(2, 0.03);
			DataFitter.fit(fnm, histGroups_p_cor_missing_n_psreg2_sec.getItem(seci), "Q");
			fnm.setLineColor(2);
			fnm.setLineWidth(3);
			fnm.setOptStat(11110);
			canpcormissingnpsreg2sec.draw(fnm, "same");
		}
		
		IndexedList<Double> meanmissingnpsreg2 = new IndexedList<>(3);
		IndexedList<Double> sigmamissingnpsreg2 = new IndexedList<>(3);
		
		IndexedList<Double> meanpcormissingnpsreg2 = new IndexedList<>(3);
		IndexedList<Double> sigmapcormissingnpsreg2 = new IndexedList<>(3);
		
		for(int seci = 0; seci < 6; seci++)
		{
			for(int theta_bini = 0; theta_bini < 9; theta_bini++)
			{
				for(int phi_bini = 0; phi_bini < 4; phi_bini++)
				{
					fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
											0.7+(histGroups_missing_n_theta_phi_bin_psreg2.getItem(seci, theta_bini, phi_bini).getMaximumBin()
													*(1.2-0.7)/100)-0.045,
											0.7+(histGroups_missing_n_theta_phi_bin_psreg2.getItem(seci, theta_bini, phi_bini).getMaximumBin()
													*(1.2-0.7)/100)+0.04);
					fnm.setParameter(0, histGroups_missing_n_theta_phi_bin_psreg2.getItem(seci, theta_bini, phi_bini).getMax()*1.2);
					fnm.setParameter(1, 0.7+(histGroups_missing_n_theta_phi_bin_psreg2.getItem(seci, theta_bini, phi_bini).getMaximumBin()
										*(1.2-0.7)/100));
					fnm.setParameter(2, 0.03);
					DataFitter.fit(fnm, histGroups_missing_n_theta_phi_bin_psreg2.getItem(seci, theta_bini, phi_bini), "Q");
					meanmissingnpsreg2.add(fnm.getParameter(1), seci, theta_bini, phi_bini);
					sigmamissingnpsreg2.add(fnm.getParameter(2), seci, theta_bini, phi_bini);
					fnm = new F1D("fnm", "[amp]*gaus(x,[mean],[sigma])",
											0.7+(histGroups_p_cor_missing_n_theta_phi_bin_psreg2.getItem(seci, theta_bini, phi_bini).getMaximumBin()
													*(1.2-0.7)/100)-0.045,
											0.7+(histGroups_p_cor_missing_n_theta_phi_bin_psreg2.getItem(seci, theta_bini, phi_bini).getMaximumBin()
													*(1.2-0.7)/100)+0.04);
					fnm.setParameter(0, histGroups_p_cor_missing_n_theta_phi_bin_psreg2.getItem(seci, theta_bini, phi_bini).getMax()*1.2);
					fnm.setParameter(1, 0.7+(histGroups_p_cor_missing_n_theta_phi_bin_psreg2.getItem(seci, theta_bini, phi_bini).getMaximumBin()
										*(1.2-0.7)/100));
					fnm.setParameter(2, 0.03);
					DataFitter.fit(fnm, histGroups_p_cor_missing_n_theta_phi_bin_psreg2.getItem(seci, theta_bini, phi_bini), "Q");
					meanpcormissingnpsreg2.add(fnm.getParameter(1), seci, theta_bini, phi_bini);
					sigmapcormissingnpsreg2.add(fnm.getParameter(2), seci, theta_bini, phi_bini);
				}
			}
		}
		
		IndexedList<GraphErrors> histerrGroups_missing_n_vs_phi_psreg2 = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> histerrGroups_p_cor_missing_n_vs_phi_psreg2 = new IndexedList<GraphErrors>(2);
		
		for(int seci = 0; seci < 6; seci++) {
			for(int theta_bini = 0; theta_bini < 9; theta_bini++) {
				GraphErrors hmissing_n_vs_phi_psreg2 = new GraphErrors();
				histerrGroups_missing_n_vs_phi_psreg2.add(hmissing_n_vs_phi_psreg2, seci, theta_bini);
				GraphErrors hp_cor_missing_n_vs_phi_psreg2 = new GraphErrors();
				histerrGroups_p_cor_missing_n_vs_phi_psreg2.add(hp_cor_missing_n_vs_phi_psreg2, seci, theta_bini);
				for(int phi_bini = 0; phi_bini < 4; phi_bini++)
				{
					if(meanmissingnpsreg2.hasItem(seci, theta_bini,  phi_bini) && sigmamissingnpsreg2.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_missing_n_vs_phi_psreg2.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]),
															meanmissingnpsreg2.getItem(seci, theta_bini,  phi_bini),
															0, sigmamissingnpsreg2.getItem(seci, theta_bini, phi_bini));
					if(meanpcormissingnpsreg2.hasItem(seci, theta_bini,  phi_bini) && sigmapcormissingnpsreg2.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_p_cor_missing_n_vs_phi_psreg2.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]),
															meanpcormissingnpsreg2.getItem(seci, theta_bini,  phi_bini),
															0, sigmapcormissingnpsreg2.getItem(seci, theta_bini, phi_bini));
				}
			}
		}
		
		for(int theta_bini = 0; theta_bini < 9; theta_bini++){
			JFrame framemissingnpsreg2phibin = new JFrame("#theta>=-10p+60: n_miss Mass vs. #phi Binned");
			framemissingnpsreg2phibin.setSize(1500, 1000);
			EmbeddedCanvas canmissingnpsreg2phibin = new EmbeddedCanvas();
			framemissingnpsreg2phibin.add(canmissingnpsreg2phibin);
			framemissingnpsreg2phibin.setLocationRelativeTo(null);
			framemissingnpsreg2phibin.setVisible(true);
			canmissingnpsreg2phibin.divide(3, 2);
			for(int seci = 0; seci < 6; seci++){
				canmissingnpsreg2phibin.cd(seci);
				canmissingnpsreg2phibin.setFont("Arial");
				histerrGroups_missing_n_vs_phi_psreg2.getItem(seci, theta_bini).setTitle("e^- Hit: " + theta_bnd[theta_bini] + "#degree < #theta <= "
																		+ theta_bnd[theta_bini+1] + " #degree, n_miss Mass vs. #phi");
				histerrGroups_missing_n_vs_phi_psreg2.getItem(seci, theta_bini).setTitleX("#phi [#degree]");
				histerrGroups_missing_n_vs_phi_psreg2.getItem(seci, theta_bini).setTitleY("Sector " + (seci+1) + " n_miss : #mu +/- #sigma [GeV]");
				histerrGroups_missing_n_vs_phi_psreg2.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_missing_n_vs_phi_psreg2.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_missing_n_vs_phi_psreg2.getItem(seci, theta_bini).setLineThickness(1);
				canmissingnpsreg2phibin.getPad(seci).setTitleFontSize(32);
				canmissingnpsreg2phibin.getPad(seci).setAxisTitleFontSize(32);
				canmissingnpsreg2phibin.getPad(seci).setAxisLabelFontSize(24);
				canmissingnpsreg2phibin.getPad(seci).setAxisRange((phi_rot[seci]-30), (phi_rot[seci]+30), 0.7, 1.2);
				histerrGroups_p_cor_missing_n_vs_phi_psreg2.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_p_cor_missing_n_vs_phi_psreg2.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_p_cor_missing_n_vs_phi_psreg2.getItem(seci, theta_bini).setMarkerColor(2);
				histerrGroups_p_cor_missing_n_vs_phi_psreg2.getItem(seci, theta_bini).setLineColor(2);
				histerrGroups_p_cor_missing_n_vs_phi_psreg2.getItem(seci, theta_bini).setLineThickness(1);			
				canmissingnpsreg2phibin.draw(histerrGroups_missing_n_vs_phi_psreg2.getItem(seci, theta_bini), "same");
				canmissingnpsreg2phibin.draw(histerrGroups_p_cor_missing_n_vs_phi_psreg2.getItem(seci, theta_bini), "same");
			}
		}
		*/
		
	// DVCS	
		
		readerdvcs.open("C:/Users/joshtanj/Documents/download/skim_6bank_epg_merged_6535eV_skim16.hipo");
		
		int dvcsEventCounter = 0;
		while(readerdvcs.hasEvent())// && dvcsEventCounter < 10000000)
		{
			dvcsEventCounter++;
			processEventdvcs(readerdvcs.getNextEvent(), dvcsEventCounter, thetaCorFac);
			if(dvcsEventCounter%500000 == 0) System.out.println("DVCS Event: " + dvcsEventCounter);
		}

		readerdvcs.close();
		
		F1D ftc;
		F1D fpm;
		F1D fXpt;
		F1D fXe;
		
		JFrame framedvcskin = new JFrame("DVCS All Sector Kinematics");
		framedvcskin.setSize(1000, 500);
		EmbeddedCanvas candvcskin = new EmbeddedCanvas();
		framedvcskin.add(candvcskin);
		framedvcskin.setLocationRelativeTo(null);
		framedvcskin.setVisible(true);
		candvcskin.divide(2, 1);
		candvcskin.cd(0);
		candvcskin.setFont("Arial");
		hdvcs_theta_vs_p_e.setTitle("e^- #theta vs. p");
		hdvcs_theta_vs_p_e.setTitleX("p [GeV]");
		hdvcs_theta_vs_p_e.setTitleY("#theta [#degree]");
		candvcskin.getPad(0).setTitleFontSize(32);
		candvcskin.getPad(0).setAxisTitleFontSize(32);
		candvcskin.getPad(0).setAxisLabelFontSize(24);
		candvcskin.draw(hdvcs_theta_vs_p_e, "same");
		candvcskin.cd(1);
		candvcskin.setFont("Arial");
		hdvcs_theta_vs_phi_e.setTitle("e^- #theta vs. #phi");
		hdvcs_theta_vs_phi_e.setTitleX("#phi [#degree]");
		hdvcs_theta_vs_phi_e.setTitleY("#theta [#degree]");
		candvcskin.getPad(1).setTitleFontSize(32);
		candvcskin.getPad(1).setAxisTitleFontSize(32);
		candvcskin.getPad(1).setAxisLabelFontSize(24);
		candvcskin.draw(hdvcs_theta_vs_phi_e, "same");
	
		JFrame framedvcsexclv = new JFrame("DVCS All Sector Exclusivity Variables");
		framedvcsexclv.setSize(1000, 1000);
		EmbeddedCanvas candvcsexclv = new EmbeddedCanvas();
		framedvcsexclv.add(candvcsexclv);
		framedvcsexclv.setLocationRelativeTo(null);
		framedvcsexclv.setVisible(true);
		candvcsexclv.divide(2, 2);
		candvcsexclv.cd(0);
		candvcsexclv.setFont("Arial");
		htheta_cone_gamma.setTitle("#Delta#theta_cone(#gamma)");
		htheta_cone_gamma.setTitleX("#Delta#theta_cone(#gamma) [#degree]");
		htheta_cone_gamma.setTitleY("Counts");
		htheta_cone_gamma.setOptStat(10);
		candvcsexclv.getPad(0).setTitleFontSize(32);
		candvcsexclv.getPad(0).setAxisTitleFontSize(32);
		candvcsexclv.getPad(0).setAxisLabelFontSize(24);
		candvcsexclv.getPad(0).setStatBoxFontSize(18);
		candvcsexclv.draw(htheta_cone_gamma, "same");
		/*
		ftc = new F1D("ftc", "[amp]*gaus(x,[mean],[sigma])", htheta_cone_gamma.getMaximumBin()*0.07, htheta_cone_gamma.getRMS()*3);
		ftc.setParameter(0, htheta_cone_gamma.getMax());
		ftc.setParLimits(0, htheta_cone_gamma.getMax()+1, htheta_cone_gamma.getMax()-1);
		ftc.setParameter(1, (htheta_cone_gamma.getMaximumBin()*0.07));
		ftc.setParLimits(1, ((htheta_cone_gamma.getMaximumBin()-1)*0.07), ((htheta_cone_gamma.getMaximumBin()+1)*0.07));
		ftc.setParameter(2, htheta_cone_gamma.getRMS()/2);
		ftc.setParLimits(2, htheta_cone_gamma.getRMS()/4, htheta_cone_gamma.getRMS());
		*/
		ftc = new F1D("ftc", "[A]*exp(-(x-[x0])/[w])", htheta_cone_gamma.getMaximumBin()*0.07, 7);
		ftc.setParameter(0, htheta_cone_gamma.getMax());
		ftc.setParLimits(0, htheta_cone_gamma.getMax()+1, htheta_cone_gamma.getMax()-1);
		ftc.setParameter(1, (htheta_cone_gamma.getMaximumBin()*0.07));
		ftc.setParLimits(1, ((htheta_cone_gamma.getMaximumBin()-1)*0.07), ((htheta_cone_gamma.getMaximumBin()+1)*0.07));
		ftc.setParameter(2, htheta_cone_gamma.getRMS()/2);
		ftc.setParLimits(2, htheta_cone_gamma.getRMS()/4, 7);
		DataFitter.fit(ftc, htheta_cone_gamma, "Q");
		ftc.setLineColor(2);
		ftc.setLineWidth(3);
		ftc.setOptStat(11110);
		candvcsexclv.draw(ftc, "same");
		candvcsexclv.cd(1);
		candvcsexclv.setFont("Arial");
		hX_eg_m.setTitle("M_X_(ep->e'#gamma)");
		hX_eg_m.setTitleX("M_X [GeV]]");
		hX_eg_m.setTitleY("Counts");
		hX_eg_m.setOptStat(10);
		candvcsexclv.getPad(1).setTitleFontSize(32);
		candvcsexclv.getPad(1).setAxisTitleFontSize(32);
		candvcsexclv.getPad(1).setAxisLabelFontSize(24);
		candvcsexclv.getPad(1).setStatBoxFontSize(18);
		candvcsexclv.draw(hX_eg_m, "same");
		fpm = new F1D("fpm", "[amp]*gaus(x,[mean],[sigma])", (hX_eg_m.getMaximumBin()*0.03-0.45), (hX_eg_m.getMaximumBin()*0.03+0.25));
		fpm.setParameter(0, hX_eg_m.getMax());
		fpm.setParameter(1, hX_eg_m.getMaximumBin()*0.03);
		fpm.setParameter(2, hX_eg_m.getRMS()/2);
		fpm.setParLimits(2, hX_eg_m.getRMS()/4, hX_eg_m.getRMS());
		DataFitter.fit(fpm, hX_eg_m, "Q");
		fpm.setLineColor(2);
		fpm.setLineWidth(3);
		fpm.setOptStat(11110);
		candvcsexclv.draw(fpm, "same");
		candvcsexclv.cd(2);
		candvcsexclv.setFont("Arial");
		hX_epg_pt.setTitle("p_#rho_X_(ep->e'p'#gamma)");
		hX_epg_pt.setTitleX("p_#rho_X [GeV]]");
		hX_epg_pt.setTitleY("Counts");
		hX_epg_pt.setOptStat(10);
		candvcsexclv.getPad(2).setTitleFontSize(32);
		candvcsexclv.getPad(2).setAxisTitleFontSize(32);
		candvcsexclv.getPad(2).setAxisLabelFontSize(24);
		candvcsexclv.getPad(2).setStatBoxFontSize(18);
		candvcsexclv.draw(hX_epg_pt, "same");
		/*
		fXpt = new F1D("fXpt", "[amp]*gaus(x,[mean],[sigma])", hX_epg_pt.getMaximumBin()*0.008, hX_epg_pt.getRMS()*3);
		fXpt.setParameter(0, hX_epg_pt.getMax());
		fXpt.setParLimits(0, hX_epg_pt.getMax()+1, hX_epg_pt.getMax()-1);
		fXpt.setParameter(1, (hX_epg_pt.getMaximumBin()*0.008));
		fXpt.setParLimits(1, ((hX_epg_pt.getMaximumBin()-1)*0.008), ((hX_epg_pt.getMaximumBin()+1)*0.008));
		fXpt.setParameter(2, hX_epg_pt.getRMS()/2);
		fXpt.setParLimits(2, hX_epg_pt.getRMS()/4, hX_epg_pt.getRMS());
		*/
		fXpt = new F1D("fXpt", "[A]*exp(-(x-[x0])/[w])", hX_epg_pt.getMaximumBin()*0.008, 0.8);
		fXpt.setParameter(0, hX_epg_pt.getMax());
		fXpt.setParLimits(0, hX_epg_pt.getMax()+1, hX_epg_pt.getMax()-1);
		fXpt.setParameter(1, (hX_epg_pt.getMaximumBin()*0.008));
		fXpt.setParLimits(1, ((hX_epg_pt.getMaximumBin()-1)*0.008), ((hX_epg_pt.getMaximumBin()+1)*0.008));
		fXpt.setParameter(2, hX_epg_pt.getRMS()/2);
		fXpt.setParLimits(2, hX_epg_pt.getRMS()/4, 0.8);
		DataFitter.fit(fXpt, hX_epg_pt, "Q");
		fXpt.setLineColor(2);
		fXpt.setLineWidth(3);
		fXpt.setOptStat(11110);
		candvcsexclv.draw(fXpt, "same");
		candvcsexclv.cd(3);
		candvcsexclv.setFont("Arial");
		hX_epg_E.setTitle("E_X_(ep->e'p'#gamma)");
		hX_epg_E.setTitleX("E_X_(ep->e'p'#gamma) [GeV]]");
		hX_epg_E.setTitleY("Counts");
		hX_epg_E.setOptStat(10);
		candvcsexclv.getPad(3).setTitleFontSize(32);
		candvcsexclv.getPad(3).setAxisTitleFontSize(32);
		candvcsexclv.getPad(3).setAxisLabelFontSize(24);
		candvcsexclv.getPad(3).setStatBoxFontSize(18);
		candvcsexclv.draw(hX_epg_E, "same");
		fXe = new F1D("fXe", "[amp]*gaus(x,[mean],[sigma])", (hX_epg_E.getMaximumBin()*0.045-2.75),
						(hX_epg_E.getMaximumBin()*0.045-2));
		fXe.setParameter(0, hX_epg_E.getMax());
		fXe.setParameter(1, (hX_epg_E.getMaximumBin()*0.045-2.25));
		fXe.setParameter(2, hX_epg_E.getRMS()/2);
		fXe.setParLimits(2, hX_epg_E.getRMS()/4, hX_epg_E.getRMS());
		DataFitter.fit(fXe, hX_epg_E, "Q");
		fXe.setLineColor(2);
		fXe.setLineWidth(3);
		fXe.setOptStat(11110);
		candvcsexclv.draw(fXe, "same");
		
		JFrame framepcordvcsexclv = new JFrame("p_e Corrected DVCS All Sector Exclusivity Variables");
		framepcordvcsexclv.setSize(1000, 1000);
		EmbeddedCanvas canpcordvcsexclv = new EmbeddedCanvas();
		framepcordvcsexclv.add(canpcordvcsexclv);
		framepcordvcsexclv.setLocationRelativeTo(null);
		framepcordvcsexclv.setVisible(true);
		canpcordvcsexclv.divide(2, 2);
		canpcordvcsexclv.cd(0);
		canpcordvcsexclv.setFont("Arial");
		hpcor_theta_cone_gamma.setTitle("#Delta#theta_cone(#gamma)");
		hpcor_theta_cone_gamma.setTitleX("#Delta#theta_cone(#gamma) [#degree]");
		hpcor_theta_cone_gamma.setTitleY("Counts");
		hpcor_theta_cone_gamma.setOptStat(10);
		canpcordvcsexclv.getPad(0).setTitleFontSize(32);
		canpcordvcsexclv.getPad(0).setAxisTitleFontSize(32);
		canpcordvcsexclv.getPad(0).setAxisLabelFontSize(24);
		canpcordvcsexclv.getPad(0).setStatBoxFontSize(18);
		canpcordvcsexclv.draw(hpcor_theta_cone_gamma, "same");
		/*
		ftc = new F1D("ftc", "[amp]*gaus(x,[mean],[sigma])", hpcor_theta_cone_gamma.getMaximumBin()*0.07,
						hpcor_theta_cone_gamma.getRMS()*3);
		ftc.setParameter(0, hpcor_theta_cone_gamma.getMax());
		ftc.setParLimits(0, hpcor_theta_cone_gamma.getMax()+1, hpcor_theta_cone_gamma.getMax()-1);
		ftc.setParameter(1, (hpcor_theta_cone_gamma.getMaximumBin()*0.07));
		ftc.setParLimits(1, ((hpcor_theta_cone_gamma.getMaximumBin()-1)*0.07), ((hpcor_theta_cone_gamma.getMaximumBin()+1)*0.07));
		ftc.setParameter(2, hpcor_theta_cone_gamma.getRMS()/2);
		ftc.setParLimits(2, hpcor_theta_cone_gamma.getRMS()/4, hpcor_theta_cone_gamma.getRMS());
		*/
		ftc = new F1D("ftc", "[A]*exp(-(x-[x0])/[w])", htheta_cone_gamma.getMaximumBin()*0.07, 7);
		ftc.setParameter(0, hpcor_theta_cone_gamma.getMax());
		ftc.setParLimits(0, hpcor_theta_cone_gamma.getMax()+1, hpcor_theta_cone_gamma.getMax()-1);
		ftc.setParameter(1, (hpcor_theta_cone_gamma.getMaximumBin()*0.07));
		ftc.setParLimits(1, ((hpcor_theta_cone_gamma.getMaximumBin()-1)*0.07), ((hpcor_theta_cone_gamma.getMaximumBin()+1)*0.07));
		ftc.setParameter(2, hpcor_theta_cone_gamma.getRMS()/2);
		ftc.setParLimits(2, hpcor_theta_cone_gamma.getRMS()/4, 7);
		DataFitter.fit(ftc, hpcor_theta_cone_gamma, "Q");
		ftc.setLineColor(2);
		ftc.setLineWidth(3);
		ftc.setOptStat(11110);
		canpcordvcsexclv.draw(ftc, "same");
		canpcordvcsexclv.cd(1);
		canpcordvcsexclv.setFont("Arial");
		hpcor_X_eg_m.setTitle("M_X_(ep->e'#gamma)");
		hpcor_X_eg_m.setTitleX("M_X [GeV]]");
		hpcor_X_eg_m.setTitleY("Counts");
		hpcor_X_eg_m.setOptStat(10);
		canpcordvcsexclv.getPad(1).setTitleFontSize(32);
		canpcordvcsexclv.getPad(1).setAxisTitleFontSize(32);
		canpcordvcsexclv.getPad(1).setAxisLabelFontSize(24);
		canpcordvcsexclv.getPad(1).setStatBoxFontSize(18);
		canpcordvcsexclv.draw(hpcor_X_eg_m, "same");
		fpm = new F1D("fpm", "[amp]*gaus(x,[mean],[sigma])", (hpcor_X_eg_m.getMaximumBin()*0.03-0.45), (hpcor_X_eg_m.getMaximumBin()*0.03+0.25));
		fpm.setParameter(0, hpcor_X_eg_m.getMax());
		fpm.setParameter(1, hpcor_X_eg_m.getMaximumBin()*0.03);
		fpm.setParameter(2, hpcor_X_eg_m.getRMS()/2);
		fpm.setParLimits(2, hpcor_X_eg_m.getRMS()/4, hpcor_X_eg_m.getRMS());
		DataFitter.fit(fpm, hpcor_X_eg_m, "Q");
		fpm.setLineColor(2);
		fpm.setLineWidth(3);
		fpm.setOptStat(11110);
		canpcordvcsexclv.draw(fpm, "same");
		canpcordvcsexclv.cd(2);
		canpcordvcsexclv.setFont("Arial");
		hpcor_X_epg_pt.setTitle("p_#rho_X_(ep->e'p'#gamma)");
		hpcor_X_epg_pt.setTitleX("p_#rho_X [GeV]]");
		hpcor_X_epg_pt.setTitleY("Counts");
		hpcor_X_epg_pt.setOptStat(10);
		canpcordvcsexclv.getPad(2).setTitleFontSize(32);
		canpcordvcsexclv.getPad(2).setAxisTitleFontSize(32);
		canpcordvcsexclv.getPad(2).setAxisLabelFontSize(24);
		canpcordvcsexclv.getPad(2).setStatBoxFontSize(18);
		canpcordvcsexclv.draw(hpcor_X_epg_pt, "same");
		/*
		fXpt = new F1D("fXpt", "[amp]*gaus(x,[mean],[sigma])", hpcor_X_epg_pt.getMaximumBin()*0.008, hpcor_X_epg_pt.getRMS()*3);
		fXpt.setParameter(0, hpcor_X_epg_pt.getMax());
		fXpt.setParLimits(0, hpcor_X_epg_pt.getMax()+1, hpcor_X_epg_pt.getMax()-1);
		fXpt.setParameter(1, (hpcor_X_epg_pt.getMaximumBin()*0.008));
		fXpt.setParLimits(1, ((hpcor_X_epg_pt.getMaximumBin()-1)*0.008), ((hpcor_X_epg_pt.getMaximumBin()+1)*0.008));
		fXpt.setParameter(2, hpcor_X_epg_pt.getRMS()/2);
		fXpt.setParLimits(2, hpcor_X_epg_pt.getRMS()/4, hpcor_X_epg_pt.getRMS());
		*/
		fXpt = new F1D("fXpt", "[A]*exp(-(x-[x0])/[w])", hpcor_X_epg_pt.getMaximumBin()*0.008, 0.8);
		fXpt.setParameter(0, hpcor_X_epg_pt.getMax());
		fXpt.setParLimits(0, hpcor_X_epg_pt.getMax()+1, hpcor_X_epg_pt.getMax()-1);
		fXpt.setParameter(1, (hpcor_X_epg_pt.getMaximumBin()*0.008));
		fXpt.setParLimits(1, ((hpcor_X_epg_pt.getMaximumBin()-1)*0.008), ((hpcor_X_epg_pt.getMaximumBin()+1)*0.008));
		fXpt.setParameter(2, hpcor_X_epg_pt.getRMS()/2);
		fXpt.setParLimits(2, hpcor_X_epg_pt.getRMS()/4, 0.8);
		DataFitter.fit(fXpt, hpcor_X_epg_pt, "Q");
		fXpt.setLineColor(2);
		fXpt.setLineWidth(3);
		fXpt.setOptStat(11110);
		canpcordvcsexclv.draw(fXpt, "same");
		canpcordvcsexclv.cd(3);
		canpcordvcsexclv.setFont("Arial");
		hpcor_X_epg_E.setTitle("E_X_(ep->e'p'#gamma)");
		hpcor_X_epg_E.setTitleX("E_X_(ep->e'p'#gamma) [GeV]]");
		hpcor_X_epg_E.setTitleY("Counts");
		hpcor_X_epg_E.setOptStat(10);
		canpcordvcsexclv.getPad(3).setTitleFontSize(32);
		canpcordvcsexclv.getPad(3).setAxisTitleFontSize(32);
		canpcordvcsexclv.getPad(3).setAxisLabelFontSize(24);
		canpcordvcsexclv.getPad(3).setStatBoxFontSize(18);
		canpcordvcsexclv.draw(hpcor_X_epg_E, "same");
		fXe = new F1D("fXe", "[amp]*gaus(x,[mean],[sigma])", (hpcor_X_epg_E.getMaximumBin()*0.045-2.75),
						(hpcor_X_epg_E.getMaximumBin()*0.045-2));
		fXe.setParameter(0, hpcor_X_epg_E.getMax());
		fXe.setParameter(1, (hpcor_X_epg_E.getMaximumBin()*0.045-2.25));
		fXe.setParameter(2, hpcor_X_epg_E.getRMS()/2);
		fXe.setParLimits(2, hpcor_X_epg_E.getRMS()/4, hpcor_X_epg_E.getRMS());
		DataFitter.fit(fXe, hpcor_X_epg_E, "Q");
		fXe.setLineColor(2);
		fXe.setLineWidth(3);
		fXe.setOptStat(11110);
		canpcordvcsexclv.draw(fXe, "same");
		
		JFrame framethetaconegammasec = new JFrame("#Delta#theta_cone(#gamma)");
		framethetaconegammasec.setSize(1500, 1000);
		EmbeddedCanvas canthetaconegammasec = new EmbeddedCanvas();
		framethetaconegammasec.add(canthetaconegammasec);
		framethetaconegammasec.setLocationRelativeTo(null);
		framethetaconegammasec.setVisible(true);
		canthetaconegammasec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canthetaconegammasec.cd(seci);
			canthetaconegammasec.setFont("Arial");
			histGroups_theta_cone_gamma_sec.getItem(seci).setTitle("Sector " + (seci+1) + " #Delta#theta_cone(#gamma)");
			histGroups_theta_cone_gamma_sec.getItem(seci).setTitleX("#Delta#theta_cone(#gamma) [#degree]");
			histGroups_theta_cone_gamma_sec.getItem(seci).setTitleY("Counts");
			histGroups_theta_cone_gamma_sec.getItem(seci).setOptStat(10);
			canthetaconegammasec.getPad(seci).setTitleFontSize(32);
			canthetaconegammasec.getPad(seci).setAxisTitleFontSize(32);
			canthetaconegammasec.getPad(seci).setAxisLabelFontSize(24);
			canthetaconegammasec.getPad(seci).setStatBoxFontSize(18);
			canthetaconegammasec.draw(histGroups_theta_cone_gamma_sec.getItem(seci), "same");
			ftc = new F1D("ftc", "[A]*exp(-(x-[x0])/[w])", htheta_cone_gamma.getMaximumBin()*0.07, 7);
			ftc.setParameter(0, histGroups_theta_cone_gamma_sec.getItem(seci).getMax());
			ftc.setParLimits(0, histGroups_theta_cone_gamma_sec.getItem(seci).getMax()+1,
									histGroups_theta_cone_gamma_sec.getItem(seci).getMax()-1);
			ftc.setParameter(1, (histGroups_theta_cone_gamma_sec.getItem(seci).getMaximumBin()*0.07));
			ftc.setParLimits(1, ((histGroups_theta_cone_gamma_sec.getItem(seci).getMaximumBin()-1)*0.07),
									((histGroups_theta_cone_gamma_sec.getItem(seci).getMaximumBin()+1)*0.07));
			ftc.setParameter(2, histGroups_theta_cone_gamma_sec.getItem(seci).getRMS()/2);
			ftc.setParLimits(2, histGroups_theta_cone_gamma_sec.getItem(seci).getRMS()/4, 7);
			DataFitter.fit(ftc, histGroups_theta_cone_gamma_sec.getItem(seci), "Q");
			ftc.setLineColor(2);
			ftc.setLineWidth(3);
			ftc.setOptStat(11110);
			canthetaconegammasec.draw(ftc, "same");
		}
		
		JFrame framepcorthetaconegammasec = new JFrame("p_e Corrected #Delta#theta_cone(#gamma)");
		framepcorthetaconegammasec.setSize(1500, 1000);
		EmbeddedCanvas canpcorthetaconegammasec = new EmbeddedCanvas();
		framepcorthetaconegammasec.add(canpcorthetaconegammasec);
		framepcorthetaconegammasec.setLocationRelativeTo(null);
		framepcorthetaconegammasec.setVisible(true);
		canpcorthetaconegammasec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canpcorthetaconegammasec.cd(seci);
			canpcorthetaconegammasec.setFont("Arial");
			histGroups_p_cor_theta_cone_gamma_sec.getItem(seci).setTitle("Sector " + (seci+1) + " #Delta#theta_cone(#gamma)");
			histGroups_p_cor_theta_cone_gamma_sec.getItem(seci).setTitleX("#Delta#theta_cone(#gamma) [#degree]");
			histGroups_p_cor_theta_cone_gamma_sec.getItem(seci).setTitleY("Counts");
			histGroups_p_cor_theta_cone_gamma_sec.getItem(seci).setOptStat(10);
			canpcorthetaconegammasec.getPad(seci).setTitleFontSize(32);
			canpcorthetaconegammasec.getPad(seci).setAxisTitleFontSize(32);
			canpcorthetaconegammasec.getPad(seci).setAxisLabelFontSize(24);
			canpcorthetaconegammasec.getPad(seci).setStatBoxFontSize(18);
			canpcorthetaconegammasec.draw(histGroups_p_cor_theta_cone_gamma_sec.getItem(seci), "same");
			ftc = new F1D("ftc", "[A]*exp(-(x-[x0])/[w])", htheta_cone_gamma.getMaximumBin()*0.07, 7);
			ftc.setParameter(0, histGroups_p_cor_theta_cone_gamma_sec.getItem(seci).getMax());
			ftc.setParLimits(0, histGroups_p_cor_theta_cone_gamma_sec.getItem(seci).getMax()+1,
									histGroups_p_cor_theta_cone_gamma_sec.getItem(seci).getMax()-1);
			ftc.setParameter(1, (histGroups_p_cor_theta_cone_gamma_sec.getItem(seci).getMaximumBin()*0.07));
			ftc.setParLimits(1, ((histGroups_p_cor_theta_cone_gamma_sec.getItem(seci).getMaximumBin()-1)*0.07),
									((histGroups_p_cor_theta_cone_gamma_sec.getItem(seci).getMaximumBin()+1)*0.07));
			ftc.setParameter(2, histGroups_p_cor_theta_cone_gamma_sec.getItem(seci).getRMS()/2);
			ftc.setParLimits(2, histGroups_p_cor_theta_cone_gamma_sec.getItem(seci).getRMS()/4, 7);
			DataFitter.fit(ftc, histGroups_p_cor_theta_cone_gamma_sec.getItem(seci), "Q");
			ftc.setLineColor(2);
			ftc.setLineWidth(3);
			ftc.setOptStat(11110);
			canpcorthetaconegammasec.draw(ftc, "same");
		}
		
		JFrame frameXegmsec = new JFrame("M_X_(ep->e'#gamma)");
		frameXegmsec.setSize(1500, 1000);
		EmbeddedCanvas canXegmsec = new EmbeddedCanvas();
		frameXegmsec.add(canXegmsec);
		frameXegmsec.setLocationRelativeTo(null);
		frameXegmsec.setVisible(true);
		canXegmsec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canXegmsec.cd(seci);
			canXegmsec.setFont("Arial");
			histGroups_X_eg_m_sec.getItem(seci).setTitle("Sector " + (seci+1) + " M_X_(ep->e'#gamma)");
			histGroups_X_eg_m_sec.getItem(seci).setTitleX("M_X [GeV]");
			histGroups_X_eg_m_sec.getItem(seci).setTitleY("Counts");
			histGroups_X_eg_m_sec.getItem(seci).setOptStat(10);
			canXegmsec.getPad(seci).setTitleFontSize(32);
			canXegmsec.getPad(seci).setAxisTitleFontSize(32);
			canXegmsec.getPad(seci).setAxisLabelFontSize(24);
			canXegmsec.getPad(seci).setStatBoxFontSize(18);
			canXegmsec.draw(histGroups_X_eg_m_sec.getItem(seci), "same");
			fpm = new F1D("fpm", "[amp]*gaus(x,[mean],[sigma])", (histGroups_X_eg_m_sec.getItem(seci).getMaximumBin()*0.03-0.45),
							(histGroups_X_eg_m_sec.getItem(seci).getMaximumBin()*0.03+0.25));
			fpm.setParameter(0, histGroups_X_eg_m_sec.getItem(seci).getMax());
			fpm.setParameter(1, histGroups_X_eg_m_sec.getItem(seci).getMaximumBin()*0.03);
			fpm.setParameter(2, histGroups_X_eg_m_sec.getItem(seci).getRMS()/2);
			fpm.setParLimits(2, histGroups_X_eg_m_sec.getItem(seci).getRMS()/4, histGroups_X_eg_m_sec.getItem(seci).getRMS());
			DataFitter.fit(fpm, histGroups_X_eg_m_sec.getItem(seci), "Q");
			fpm.setLineColor(2);
			fpm.setLineWidth(3);
			fpm.setOptStat(11110);
			canXegmsec.draw(fpm, "same");
		}
		
		JFrame framepcorXegmsec = new JFrame("p_e Corrected M_X_(ep->e'#gamma)");
		framepcorXegmsec.setSize(1500, 1000);
		EmbeddedCanvas canpcorXegmsec = new EmbeddedCanvas();
		framepcorXegmsec.add(canpcorXegmsec);
		framepcorXegmsec.setLocationRelativeTo(null);
		framepcorXegmsec.setVisible(true);
		canpcorXegmsec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canpcorXegmsec.cd(seci);
			canpcorXegmsec.setFont("Arial");
			histGroups_p_cor_X_eg_m_sec.getItem(seci).setTitle("Sector " + (seci+1) + " M_X_(ep->e'#gamma)");
			histGroups_p_cor_X_eg_m_sec.getItem(seci).setTitleX("M_X [GeV]");
			histGroups_p_cor_X_eg_m_sec.getItem(seci).setTitleY("Counts");
			histGroups_p_cor_X_eg_m_sec.getItem(seci).setOptStat(10);
			canpcorXegmsec.getPad(seci).setTitleFontSize(32);
			canpcorXegmsec.getPad(seci).setAxisTitleFontSize(32);
			canpcorXegmsec.getPad(seci).setAxisLabelFontSize(24);
			canpcorXegmsec.getPad(seci).setStatBoxFontSize(18);
			canpcorXegmsec.draw(histGroups_p_cor_X_eg_m_sec.getItem(seci), "same");
			fpm = new F1D("fpm", "[amp]*gaus(x,[mean],[sigma])", (histGroups_p_cor_X_eg_m_sec.getItem(seci).getMaximumBin()*0.03-0.45),
							(histGroups_p_cor_X_eg_m_sec.getItem(seci).getMaximumBin()*0.03+0.25));
			fpm.setParameter(0, histGroups_p_cor_X_eg_m_sec.getItem(seci).getMax());
			fpm.setParameter(1, histGroups_p_cor_X_eg_m_sec.getItem(seci).getMaximumBin()*0.03);
			fpm.setParameter(2, histGroups_p_cor_X_eg_m_sec.getItem(seci).getRMS()/2);
			fpm.setParLimits(2, histGroups_p_cor_X_eg_m_sec.getItem(seci).getRMS()/4, histGroups_p_cor_X_eg_m_sec.getItem(seci).getRMS());
			DataFitter.fit(fpm, histGroups_p_cor_X_eg_m_sec.getItem(seci), "Q");
			fpm.setLineColor(2);
			fpm.setLineWidth(3);
			fpm.setOptStat(11110);
			canpcorXegmsec.draw(fpm, "same");
		}
		
		JFrame frameXepgptsec = new JFrame(" p_#rho_X_(ep->e'p'#gamma)");
		frameXepgptsec.setSize(1500, 1000);
		EmbeddedCanvas canXepgptsec = new EmbeddedCanvas();
		frameXepgptsec.add(canXepgptsec);
		frameXepgptsec.setLocationRelativeTo(null);
		frameXepgptsec.setVisible(true);
		canXepgptsec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canXepgptsec.cd(seci);
			canXepgptsec.setFont("Arial");
			histGroups_X_epg_pt_sec.getItem(seci).setTitle("Sector " + (seci+1) + " p_#rho_X_(ep->e'p'#gamma)");
			histGroups_X_epg_pt_sec.getItem(seci).setTitleX("M_X [GeV]");
			histGroups_X_epg_pt_sec.getItem(seci).setTitleY("Counts");
			histGroups_X_epg_pt_sec.getItem(seci).setOptStat(10);
			canXepgptsec.getPad(seci).setTitleFontSize(32);
			canXepgptsec.getPad(seci).setAxisTitleFontSize(32);
			canXepgptsec.getPad(seci).setAxisLabelFontSize(24);
			canXepgptsec.getPad(seci).setStatBoxFontSize(18);
			canXepgptsec.draw(histGroups_X_epg_pt_sec.getItem(seci), "same");
			fXpt = new F1D("fXpt", "[A]*exp(-(x-[x0])/[w])",  histGroups_X_epg_pt_sec.getItem(seci).getMaximumBin()*0.008, 0.8);
			fXpt.setParameter(0,  histGroups_X_epg_pt_sec.getItem(seci).getMax());
			fXpt.setParLimits(0,  histGroups_X_epg_pt_sec.getItem(seci).getMax()+1,  histGroups_X_epg_pt_sec.getItem(seci).getMax()-1);
			fXpt.setParameter(1, ( histGroups_X_epg_pt_sec.getItem(seci).getMaximumBin()*0.008));
			fXpt.setParLimits(1, (( histGroups_X_epg_pt_sec.getItem(seci).getMaximumBin()-1)*0.008), (( histGroups_X_epg_pt_sec.getItem(seci).getMaximumBin()+1)*0.008));
			fXpt.setParameter(2,  histGroups_X_epg_pt_sec.getItem(seci).getRMS()/2);
			fXpt.setParLimits(2,  histGroups_X_epg_pt_sec.getItem(seci).getRMS()/4, 0.8);
			DataFitter.fit(fXpt, histGroups_X_epg_pt_sec.getItem(seci), "Q");
			fXpt.setLineColor(2);
			fXpt.setLineWidth(3);
			fXpt.setOptStat(11110);
			canXepgptsec.draw(fXpt, "same");
		}
		
		JFrame framepcorXepgptsec = new JFrame("p_e Corrected p_#rho_X_(ep->e'p'#gamma)");
		framepcorXepgptsec.setSize(1500, 1000);
		EmbeddedCanvas canpcorXepgptsec = new EmbeddedCanvas();
		framepcorXepgptsec.add(canpcorXepgptsec);
		framepcorXepgptsec.setLocationRelativeTo(null);
		framepcorXepgptsec.setVisible(true);
		canpcorXepgptsec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canpcorXepgptsec.cd(seci);
			canpcorXepgptsec.setFont("Arial");
			histGroups_p_cor_X_epg_pt_sec.getItem(seci).setTitle("Sector " + (seci+1) + " p_#rho_X_(ep->e'p'#gamma)");
			histGroups_p_cor_X_epg_pt_sec.getItem(seci).setTitleX("M_X [GeV]");
			histGroups_p_cor_X_epg_pt_sec.getItem(seci).setTitleY("Counts");
			histGroups_p_cor_X_epg_pt_sec.getItem(seci).setOptStat(10);
			canpcorXepgptsec.getPad(seci).setTitleFontSize(32);
			canpcorXepgptsec.getPad(seci).setAxisTitleFontSize(32);
			canpcorXepgptsec.getPad(seci).setAxisLabelFontSize(24);
			canpcorXepgptsec.getPad(seci).setStatBoxFontSize(18);
			canpcorXepgptsec.draw(histGroups_p_cor_X_epg_pt_sec.getItem(seci), "same");
			fXpt = new F1D("fXpt", "[A]*exp(-(x-[x0])/[w])",  histGroups_p_cor_X_epg_pt_sec.getItem(seci).getMaximumBin()*0.008, 0.8);
			fXpt.setParameter(0,  histGroups_p_cor_X_epg_pt_sec.getItem(seci).getMax());
			fXpt.setParLimits(0,  histGroups_p_cor_X_epg_pt_sec.getItem(seci).getMax()+1,  histGroups_p_cor_X_epg_pt_sec.getItem(seci).getMax()-1);
			fXpt.setParameter(1, ( histGroups_p_cor_X_epg_pt_sec.getItem(seci).getMaximumBin()*0.008));
			fXpt.setParLimits(1, (( histGroups_p_cor_X_epg_pt_sec.getItem(seci).getMaximumBin()-1)*0.008),
								(( histGroups_p_cor_X_epg_pt_sec.getItem(seci).getMaximumBin()+1)*0.008));
			fXpt.setParameter(2,  histGroups_p_cor_X_epg_pt_sec.getItem(seci).getRMS()/2);
			fXpt.setParLimits(2,  histGroups_p_cor_X_epg_pt_sec.getItem(seci).getRMS()/4, 0.8);
			DataFitter.fit(fXpt, histGroups_p_cor_X_epg_pt_sec.getItem(seci), "Q");
			fXpt.setLineColor(2);
			fXpt.setLineWidth(3);
			fXpt.setOptStat(11110);
			canpcorXepgptsec.draw(fXpt, "same");
		}
		
		JFrame frameXepgEsec = new JFrame("E_X_(ep->e'p'#gamma)");
		frameXepgEsec.setSize(1500, 1000);
		EmbeddedCanvas canXepgEsec = new EmbeddedCanvas();
		frameXepgEsec.add(canXepgEsec);
		frameXepgEsec.setLocationRelativeTo(null);
		frameXepgEsec.setVisible(true);
		canXepgEsec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canXepgEsec.cd(seci);
			canXepgEsec.setFont("Arial");
			histGroups_X_epg_E_sec.getItem(seci).setTitle("Sector " + (seci+1) + " E_X_(ep->e'p'#gamma)");
			histGroups_X_epg_E_sec.getItem(seci).setTitleX("E_X_(ep->e'p'#gamma) [GeV]");
			histGroups_X_epg_E_sec.getItem(seci).setTitleY("Counts");
			histGroups_X_epg_E_sec.getItem(seci).setOptStat(10);
			canXepgEsec.getPad(seci).setTitleFontSize(32);
			canXepgEsec.getPad(seci).setAxisTitleFontSize(32);
			canXepgEsec.getPad(seci).setAxisLabelFontSize(24);
			canXepgEsec.getPad(seci).setStatBoxFontSize(18);
			canXepgEsec.draw(histGroups_X_epg_E_sec.getItem(seci), "same");
			fXe = new F1D("fXe", "[amp]*gaus(x,[mean],[sigma])", (histGroups_X_epg_E_sec.getItem(seci).getMaximumBin()*0.045-2.75),
							(histGroups_X_epg_E_sec.getItem(seci).getMaximumBin()*0.045-2));
			fXe.setParameter(0, histGroups_X_epg_E_sec.getItem(seci).getMax());
			fXe.setParameter(1, (histGroups_X_epg_E_sec.getItem(seci).getMaximumBin()*0.045-2.25));
			fXe.setParameter(2, histGroups_X_epg_E_sec.getItem(seci).getRMS()/2);
			fXe.setParLimits(2, histGroups_X_epg_E_sec.getItem(seci).getRMS()/4, histGroups_X_epg_E_sec.getItem(seci).getRMS());
			DataFitter.fit(fXe, histGroups_X_epg_E_sec.getItem(seci), "Q");
			fXe.setLineColor(2);
			fXe.setLineWidth(3);
			fXe.setOptStat(11110);
			canXepgEsec.draw(fXe, "same");
		}
		
		JFrame framepcorXepgEsec = new JFrame("p_e Corrected E_X_(ep->e'p'#gamma)");
		framepcorXepgEsec.setSize(1500, 1000);
		EmbeddedCanvas canpcorXepgEsec = new EmbeddedCanvas();
		framepcorXepgEsec.add(canpcorXepgEsec);
		framepcorXepgEsec.setLocationRelativeTo(null);
		framepcorXepgEsec.setVisible(true);
		canpcorXepgEsec.divide(3, 2);
		for(int seci = 0; seci < 6; seci++){
			canpcorXepgEsec.cd(seci);
			canpcorXepgEsec.setFont("Arial");
			histGroups_p_cor_X_epg_E_sec.getItem(seci).setTitle("Sector " + (seci+1) + " E_X_(ep->e'p'#gamma)");
			histGroups_p_cor_X_epg_E_sec.getItem(seci).setTitleX("E_X_(ep->e'p'#gamma) [GeV]");
			histGroups_p_cor_X_epg_E_sec.getItem(seci).setTitleY("Counts");
			histGroups_p_cor_X_epg_E_sec.getItem(seci).setOptStat(10);
			canpcorXepgEsec.getPad(seci).setTitleFontSize(32);
			canpcorXepgEsec.getPad(seci).setAxisTitleFontSize(32);
			canpcorXepgEsec.getPad(seci).setAxisLabelFontSize(24);
			canpcorXepgEsec.getPad(seci).setStatBoxFontSize(18);
			canpcorXepgEsec.draw(histGroups_p_cor_X_epg_E_sec.getItem(seci), "same");
			fXe = new F1D("fXe", "[amp]*gaus(x,[mean],[sigma])", (histGroups_p_cor_X_epg_E_sec.getItem(seci).getMaximumBin()*0.045-2.75),
							(histGroups_p_cor_X_epg_E_sec.getItem(seci).getMaximumBin()*0.045-2));
			fXe.setParameter(0, histGroups_p_cor_X_epg_E_sec.getItem(seci).getMax());
			fXe.setParameter(1, (histGroups_p_cor_X_epg_E_sec.getItem(seci).getMaximumBin()*0.045-2.25));
			fXe.setParameter(2, histGroups_p_cor_X_epg_E_sec.getItem(seci).getRMS()/2);
			fXe.setParLimits(2, histGroups_p_cor_X_epg_E_sec.getItem(seci).getRMS()/4, histGroups_p_cor_X_epg_E_sec.getItem(seci).getRMS());
			DataFitter.fit(fXe, histGroups_p_cor_X_epg_E_sec.getItem(seci), "Q");
			fXe.setLineColor(2);
			fXe.setLineWidth(3);
			fXe.setOptStat(11110);
			canpcorXepgEsec.draw(fXe, "same");
		}
		
		IndexedList<Double> meanthetaconegamma = new IndexedList<>(3);
		IndexedList<Double> sigmathetaconegamma = new IndexedList<>(3);
		
		IndexedList<Double> meanpcorthetaconegamma = new IndexedList<>(3);
		IndexedList<Double> sigmapcorthetaconegamma = new IndexedList<>(3);
		
		IndexedList<Double> meanXegm = new IndexedList<>(3);
		IndexedList<Double> sigmaXegm = new IndexedList<>(3);
		
		IndexedList<Double> meanpcorXegm = new IndexedList<>(3);
		IndexedList<Double> sigmapcorXegm = new IndexedList<>(3);
		
		IndexedList<Double> meanXepgpt = new IndexedList<>(3);
		IndexedList<Double> sigmaXepgpt = new IndexedList<>(3);
		
		IndexedList<Double> meanpcorXepgpt = new IndexedList<>(3);
		IndexedList<Double> sigmapcorXepgpt = new IndexedList<>(3);
		
		IndexedList<Double> meanXepgE = new IndexedList<>(3);
		IndexedList<Double> sigmaXepgE = new IndexedList<>(3);
		
		IndexedList<Double> meanpcorXepgE = new IndexedList<>(3);
		IndexedList<Double> sigmapcorXepgE = new IndexedList<>(3);
		
		for(int seci = 0; seci < 6; seci++)
		{
			for(int theta_bini = 0; theta_bini < 9; theta_bini++)
			{
				for(int phi_bini = 0; phi_bini < 4; phi_bini++)
				{
					ftc = new F1D("ftc", "[A]*exp(-(x-[x0])/[w])",
									histGroups_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.07, 7);
					ftc.setParameter(0, histGroups_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax());
					ftc.setParLimits(0, histGroups_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax()+1,
										histGroups_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax()-1);
					ftc.setParameter(1, (histGroups_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.07));
					ftc.setParLimits(1, ((histGroups_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()-1)*0.07),
										((histGroups_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()+1)*0.07));
					ftc.setParameter(2, histGroups_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/2);
					ftc.setParLimits(2, histGroups_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/4, 7);
					DataFitter.fit(ftc, histGroups_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini), "Q");
					meanthetaconegamma.add(ftc.getParameter(1), seci, theta_bini, phi_bini);
					sigmathetaconegamma.add(ftc.getParameter(2), seci, theta_bini, phi_bini);
					ftc = new F1D("ftc", "[A]*exp(-(x-[x0])/[w])",
									histGroups_p_cor_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.07, 7);
					ftc.setParameter(0, histGroups_p_cor_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax());
					ftc.setParLimits(0, histGroups_p_cor_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax()+1,
										histGroups_p_cor_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax()-1);
					ftc.setParameter(1, (histGroups_p_cor_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.07));
					ftc.setParLimits(1, ((histGroups_p_cor_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()-1)
											*0.07),
										((histGroups_p_cor_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()+1)
											*0.07));
					ftc.setParameter(2, histGroups_p_cor_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/2);
					ftc.setParLimits(2, histGroups_p_cor_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/4, 7);
					DataFitter.fit(ftc, histGroups_p_cor_theta_cone_gamma_theta_phi_bin.getItem(seci, theta_bini, phi_bini), "Q");
					meanpcorthetaconegamma.add(ftc.getParameter(1), seci, theta_bini, phi_bini);
					sigmapcorthetaconegamma.add(ftc.getParameter(2), seci, theta_bini, phi_bini);
					
					fpm = new F1D("fpm", "[amp]*gaus(x,[mean],[sigma])",
									(histGroups_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.03-0.45),
									(histGroups_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.03+0.25));
					fpm.setParameter(0, histGroups_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax());
					fpm.setParameter(1, histGroups_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.03);
					fpm.setParameter(2, histGroups_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/2);
					fpm.setParLimits(2, histGroups_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/4,
										histGroups_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS());
					DataFitter.fit(fpm, histGroups_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini), "Q");
					meanXegm.add(fpm.getParameter(1), seci, theta_bini, phi_bini);
					sigmaXegm.add(fpm.getParameter(2), seci, theta_bini, phi_bini);
					fpm = new F1D("fpm", "[amp]*gaus(x,[mean],[sigma])",
									(histGroups_p_cor_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.03-0.45),
									(histGroups_p_cor_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.03+0.25));
					fpm.setParameter(0, histGroups_p_cor_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax());
					fpm.setParameter(1, histGroups_p_cor_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.03);
					fpm.setParameter(2, histGroups_p_cor_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/2);
					fpm.setParLimits(2, histGroups_p_cor_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/4,
										histGroups_p_cor_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS());
					DataFitter.fit(fpm, histGroups_p_cor_X_eg_m_theta_phi_bin.getItem(seci, theta_bini, phi_bini), "Q");
					meanpcorXegm.add(fpm.getParameter(1), seci, theta_bini, phi_bini);
					sigmapcorXegm.add(fpm.getParameter(2), seci, theta_bini, phi_bini);
					
					fXpt = new F1D("fXpt", "[A]*exp(-(x-[x0])/[w])",
									histGroups_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.008, 0.8);
					fXpt.setParameter(0, histGroups_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax());
					fXpt.setParLimits(0, histGroups_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax()+1,
										histGroups_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax()-1);
					fXpt.setParameter(1, (histGroups_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.008));
					fXpt.setParLimits(1, ((histGroups_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()-1)*0.008),
										((histGroups_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()+1)*0.008));
					fXpt.setParameter(2, histGroups_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/2);
					fXpt.setParLimits(2, histGroups_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/4, 0.8);
					DataFitter.fit(fXpt, histGroups_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini), "Q");
					meanXepgpt.add(fXpt.getParameter(1), seci, theta_bini, phi_bini);
					sigmaXepgpt.add(fXpt.getParameter(2), seci, theta_bini, phi_bini);
					fXpt = new F1D("fXpt", "[A]*exp(-(x-[x0])/[w])",
									histGroups_p_cor_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.008, 0.8);
					fXpt.setParameter(0, histGroups_p_cor_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax());
					fXpt.setParLimits(0, histGroups_p_cor_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax()+1,
										histGroups_p_cor_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax()-1);
					fXpt.setParameter(1, (histGroups_p_cor_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.008));
					fXpt.setParLimits(1, ((histGroups_p_cor_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()-1)*0.008),
										((histGroups_p_cor_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()+1)*0.008));
					fXpt.setParameter(2, histGroups_p_cor_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/2);
					fXpt.setParLimits(2, histGroups_p_cor_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/4, 0.8);
					DataFitter.fit(fXpt, histGroups_p_cor_X_epg_pt_theta_phi_bin.getItem(seci, theta_bini, phi_bini), "Q");
					meanpcorXepgpt.add(fXpt.getParameter(1), seci, theta_bini, phi_bini);
					sigmapcorXepgpt.add(fXpt.getParameter(2), seci, theta_bini, phi_bini);
					
					fXe = new F1D("fXe", "[amp]*gaus(x,[mean],[sigma])",
									(histGroups_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.045-2.75),
									(histGroups_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.045-2));
					fXe.setParameter(0, histGroups_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax());
					fXe.setParameter(1, (histGroups_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.045-2.25));
					fXe.setParameter(2, histGroups_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/2);
					fXe.setParLimits(2, histGroups_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/4,
										histGroups_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS());
					DataFitter.fit(fXe, histGroups_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini), "Q");
					meanXepgE.add(fXe.getParameter(1), seci, theta_bini, phi_bini);
					sigmaXepgE.add(fXe.getParameter(2), seci, theta_bini, phi_bini);
					fXe = new F1D("fXe", "[amp]*gaus(x,[mean],[sigma])",
									(histGroups_p_cor_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.045-2.75),
									(histGroups_p_cor_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.045-2));
					fXe.setParameter(0, histGroups_p_cor_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMax());
					fXe.setParameter(1, (histGroups_p_cor_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getMaximumBin()*0.045-2.25));
					fXe.setParameter(2, histGroups_p_cor_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/2);
					fXe.setParLimits(2, histGroups_p_cor_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS()/4,
										histGroups_p_cor_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini).getRMS());
					DataFitter.fit(fXe, histGroups_p_cor_X_epg_E_theta_phi_bin.getItem(seci, theta_bini, phi_bini), "Q");
					meanpcorXepgE.add(fXe.getParameter(1), seci, theta_bini, phi_bini);
					sigmapcorXepgE.add(fXe.getParameter(2), seci, theta_bini, phi_bini);
				}
			}
		}
		
		IndexedList<GraphErrors> histerrGroups_theta_cone_gamma_vs_phi = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> histerrGroups_p_cor_theta_cone_gamma_vs_phi = new IndexedList<GraphErrors>(2);
		
		IndexedList<GraphErrors> histerrGroups_X_eg_m_vs_phi = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> histerrGroups_p_cor_X_eg_m_vs_phi = new IndexedList<GraphErrors>(2);
		
		IndexedList<GraphErrors> histerrGroups_X_epg_pt_vs_phi = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> histerrGroups_p_cor_X_epg_pt_vs_phi = new IndexedList<GraphErrors>(2);
		
		IndexedList<GraphErrors> histerrGroups_X_epg_E_vs_phi = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> histerrGroups_p_cor_X_epg_E_vs_phi = new IndexedList<GraphErrors>(2);
		
		for(int seci = 0; seci < 6; seci++) {
			for(int theta_bini = 0; theta_bini < 9; theta_bini++) {
				GraphErrors htheta_cone_gamma_vs_phi = new GraphErrors();
				histerrGroups_theta_cone_gamma_vs_phi.add(htheta_cone_gamma_vs_phi, seci, theta_bini);
				GraphErrors hp_cor_theta_cone_gamma_vs_phi = new GraphErrors();
				histerrGroups_p_cor_theta_cone_gamma_vs_phi.add(hp_cor_theta_cone_gamma_vs_phi, seci, theta_bini);
				
				GraphErrors hX_eg_m_vs_phi = new GraphErrors();
				histerrGroups_X_eg_m_vs_phi.add(hX_eg_m_vs_phi, seci, theta_bini);
				GraphErrors hp_cor_X_eg_m_vs_phi = new GraphErrors();
				histerrGroups_p_cor_X_eg_m_vs_phi.add(hp_cor_X_eg_m_vs_phi, seci, theta_bini);
				
				GraphErrors hX_epg_pt_vs_phi = new GraphErrors();
				histerrGroups_X_epg_pt_vs_phi.add(hX_epg_pt_vs_phi, seci, theta_bini);
				GraphErrors hp_cor_X_epg_pt_vs_phi = new GraphErrors();
				histerrGroups_p_cor_X_epg_pt_vs_phi.add(hp_cor_X_epg_pt_vs_phi, seci, theta_bini);
				
				GraphErrors hX_epg_E_vs_phi = new GraphErrors();
				histerrGroups_X_epg_E_vs_phi.add(hX_epg_E_vs_phi, seci, theta_bini);
				GraphErrors hp_cor_X_epg_E_vs_phi = new GraphErrors();
				histerrGroups_p_cor_X_epg_E_vs_phi.add(hp_cor_X_epg_E_vs_phi, seci, theta_bini);
				
				for(int phi_bini = 0; phi_bini < 4; phi_bini++)
				{
					if(meanthetaconegamma.hasItem(seci, theta_bini,  phi_bini) && sigmathetaconegamma.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_theta_cone_gamma_vs_phi.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]),
															meanthetaconegamma.getItem(seci, theta_bini,  phi_bini),
															0, sigmathetaconegamma.getItem(seci, theta_bini, phi_bini));
					if(meanpcorthetaconegamma.hasItem(seci, theta_bini,  phi_bini) && sigmapcorthetaconegamma.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_p_cor_theta_cone_gamma_vs_phi.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]),
															meanpcorthetaconegamma.getItem(seci, theta_bini,  phi_bini),
															0, sigmapcorthetaconegamma.getItem(seci, theta_bini, phi_bini));
					
					if(meanXegm.hasItem(seci, theta_bini,  phi_bini) && sigmaXegm.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_X_eg_m_vs_phi.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]),
															meanXegm.getItem(seci, theta_bini,  phi_bini),
															0, sigmaXegm.getItem(seci, theta_bini, phi_bini));
					if(meanpcorXegm.hasItem(seci, theta_bini,  phi_bini) && sigmapcorXegm.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_p_cor_X_eg_m_vs_phi.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]),
															meanpcorXegm.getItem(seci, theta_bini,  phi_bini),
															0, sigmapcorXegm.getItem(seci, theta_bini, phi_bini));
					
					if(meanXepgpt.hasItem(seci, theta_bini,  phi_bini) && sigmaXepgpt.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_X_epg_pt_vs_phi.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]),
															meanXepgpt.getItem(seci, theta_bini,  phi_bini),
															0, sigmaXepgpt.getItem(seci, theta_bini, phi_bini));
					if(meanpcorXepgpt.hasItem(seci, theta_bini,  phi_bini) && sigmapcorXepgpt.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_p_cor_X_epg_pt_vs_phi.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]),
															meanpcorXepgpt.getItem(seci, theta_bini,  phi_bini),
															0, sigmapcorXepgpt.getItem(seci, theta_bini, phi_bini));
					
					if(meanXepgE.hasItem(seci, theta_bini,  phi_bini) && sigmaXepgE.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_X_epg_E_vs_phi.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]),
															meanXepgE.getItem(seci, theta_bini,  phi_bini),
															0, sigmaXepgE.getItem(seci, theta_bini, phi_bini));
					if(meanpcorXepgE.hasItem(seci, theta_bini,  phi_bini) && sigmapcorXepgE.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_p_cor_X_epg_E_vs_phi.getItem(seci, theta_bini)
												.addPoint((phirotba[phi_bini]+phi_rot[seci]),
															meanpcorXepgE.getItem(seci, theta_bini,  phi_bini),
															0, sigmapcorXepgE.getItem(seci, theta_bini, phi_bini));
				}
			}
		}
		
		for(int theta_bini = 0; theta_bini < 9; theta_bini++){
			JFrame framethetaconegammaphibin = new JFrame("#Delta#theta_cone(#gamma) vs. #phi Binned");
			framethetaconegammaphibin.setSize(1500, 1000);
			EmbeddedCanvas canthetaconegammaphibin = new EmbeddedCanvas();
			framethetaconegammaphibin.add(canthetaconegammaphibin);
			framethetaconegammaphibin.setLocationRelativeTo(null);
			framethetaconegammaphibin.setVisible(true);
			canthetaconegammaphibin.divide(3, 2);
			for(int seci = 0; seci < 6; seci++){
				canthetaconegammaphibin.cd(seci);
				canthetaconegammaphibin.setFont("Arial");
				histerrGroups_theta_cone_gamma_vs_phi.getItem(seci, theta_bini).setTitle("e^- Hit: " + theta_bnd[theta_bini] + "#degree < #theta <= "
																		+ theta_bnd[theta_bini+1] + " #degree, #Delta#theta_cone(#gamma) vs. #phi");
				histerrGroups_theta_cone_gamma_vs_phi.getItem(seci, theta_bini).setTitleX("#phi [#degree]");
				histerrGroups_theta_cone_gamma_vs_phi.getItem(seci, theta_bini).setTitleY("Sector " + (seci+1)
																							+ " #Delta#theta_cone(#gamma) : x_0 +/- #w [#degree]");
				histerrGroups_theta_cone_gamma_vs_phi.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_theta_cone_gamma_vs_phi.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_theta_cone_gamma_vs_phi.getItem(seci, theta_bini).setLineThickness(1);
				canthetaconegammaphibin.getPad(seci).setTitleFontSize(32);
				canthetaconegammaphibin.getPad(seci).setAxisTitleFontSize(32);
				canthetaconegammaphibin.getPad(seci).setAxisLabelFontSize(24);
				canthetaconegammaphibin.getPad(seci).setAxisRange((phi_rot[seci]-30), (phi_rot[seci]+30), 0, 7);
				histerrGroups_p_cor_theta_cone_gamma_vs_phi.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_p_cor_theta_cone_gamma_vs_phi.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_p_cor_theta_cone_gamma_vs_phi.getItem(seci, theta_bini).setMarkerColor(2);
				histerrGroups_p_cor_theta_cone_gamma_vs_phi.getItem(seci, theta_bini).setLineColor(2);
				histerrGroups_p_cor_theta_cone_gamma_vs_phi.getItem(seci, theta_bini).setLineThickness(1);
				canthetaconegammaphibin.draw(histerrGroups_theta_cone_gamma_vs_phi.getItem(seci, theta_bini), "same");
				canthetaconegammaphibin.draw(histerrGroups_p_cor_theta_cone_gamma_vs_phi.getItem(seci, theta_bini), "same");
			}
		}
		
		for(int theta_bini = 0; theta_bini < 9; theta_bini++){
			JFrame frameXegmphibin = new JFrame("M_X_(ep->e'#gamma) vs. #phi Binned");
			frameXegmphibin.setSize(1500, 1000);
			EmbeddedCanvas canXegmphibin = new EmbeddedCanvas();
			frameXegmphibin.add(canXegmphibin);
			frameXegmphibin.setLocationRelativeTo(null);
			frameXegmphibin.setVisible(true);
			canXegmphibin.divide(3, 2);
			for(int seci = 0; seci < 6; seci++){
				canXegmphibin.cd(seci);
				canXegmphibin.setFont("Arial");
				histerrGroups_X_eg_m_vs_phi.getItem(seci, theta_bini).setTitle("e^- Hit: " + theta_bnd[theta_bini] + "#degree < #theta <= "
																		+ theta_bnd[theta_bini+1] + " #degree, M_X_(ep->e'#gamma) vs. #phi");
				histerrGroups_X_eg_m_vs_phi.getItem(seci, theta_bini).setTitleX("#phi [#degree]");
				histerrGroups_X_eg_m_vs_phi.getItem(seci, theta_bini).setTitleY("Sector " + (seci+1) + " M_X : #mu +/- #sigma [GeV]");
				histerrGroups_X_eg_m_vs_phi.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_X_eg_m_vs_phi.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_X_eg_m_vs_phi.getItem(seci, theta_bini).setLineThickness(1);
				canXegmphibin.getPad(seci).setTitleFontSize(32);
				canXegmphibin.getPad(seci).setAxisTitleFontSize(32);
				canXegmphibin.getPad(seci).setAxisLabelFontSize(24);
				canXegmphibin.getPad(seci).setAxisRange((phi_rot[seci]-30), (phi_rot[seci]+30), 0, 3);
				histerrGroups_p_cor_X_eg_m_vs_phi.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_p_cor_X_eg_m_vs_phi.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_p_cor_X_eg_m_vs_phi.getItem(seci, theta_bini).setMarkerColor(2);
				histerrGroups_p_cor_X_eg_m_vs_phi.getItem(seci, theta_bini).setLineColor(2);
				histerrGroups_p_cor_X_eg_m_vs_phi.getItem(seci, theta_bini).setLineThickness(1);
				canXegmphibin.draw(histerrGroups_X_eg_m_vs_phi.getItem(seci, theta_bini), "same");
				canXegmphibin.draw(histerrGroups_p_cor_X_eg_m_vs_phi.getItem(seci, theta_bini), "same");
			}
		}
		
		for(int theta_bini = 0; theta_bini < 9; theta_bini++){
			JFrame frameXepgptphibin = new JFrame("p_#rho_X_(ep->e'p'#gamma) vs. #phi Binned");
			frameXepgptphibin.setSize(1500, 1000);
			EmbeddedCanvas canXepgptphibin = new EmbeddedCanvas();
			frameXepgptphibin.add(canXepgptphibin);
			frameXepgptphibin.setLocationRelativeTo(null);
			frameXepgptphibin.setVisible(true);
			canXepgptphibin.divide(3, 2);
			for(int seci = 0; seci < 6; seci++){
				canXepgptphibin.cd(seci);
				canXepgptphibin.setFont("Arial");
				histerrGroups_X_epg_pt_vs_phi.getItem(seci, theta_bini).setTitle("e^- Hit: " + theta_bnd[theta_bini] + "#degree < #theta <= "
																		+ theta_bnd[theta_bini+1] + " p_#rho_X_(ep->e'p'#gamma), n_miss Mass vs. #phi");
				histerrGroups_X_epg_pt_vs_phi.getItem(seci, theta_bini).setTitleX("#phi [#degree]");
				histerrGroups_X_epg_pt_vs_phi.getItem(seci, theta_bini).setTitleY("Sector " + (seci+1) + "  p_#rho_X : x_0 +/- w [GeV]");
				histerrGroups_X_epg_pt_vs_phi.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_X_epg_pt_vs_phi.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_X_epg_pt_vs_phi.getItem(seci, theta_bini).setLineThickness(1);
				canXepgptphibin.getPad(seci).setTitleFontSize(32);
				canXepgptphibin.getPad(seci).setAxisTitleFontSize(32);
				canXepgptphibin.getPad(seci).setAxisLabelFontSize(24);
				canXepgptphibin.getPad(seci).setAxisRange((phi_rot[seci]-30), (phi_rot[seci]+30), 0, 0.8);
				histerrGroups_p_cor_X_epg_pt_vs_phi.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_p_cor_X_epg_pt_vs_phi.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_p_cor_X_epg_pt_vs_phi.getItem(seci, theta_bini).setMarkerColor(2);
				histerrGroups_p_cor_X_epg_pt_vs_phi.getItem(seci, theta_bini).setLineColor(2);
				histerrGroups_p_cor_X_epg_pt_vs_phi.getItem(seci, theta_bini).setLineThickness(1);
				canXepgptphibin.draw(histerrGroups_X_epg_pt_vs_phi.getItem(seci, theta_bini), "same");
				canXepgptphibin.draw(histerrGroups_p_cor_X_epg_pt_vs_phi.getItem(seci, theta_bini), "same");
			}
		}
		
		for(int theta_bini = 0; theta_bini < 9; theta_bini++){
			JFrame frameXepgEphibin = new JFrame("E_X_(ep->e'p'#gamma) vs. #phi Binned");
			frameXepgEphibin.setSize(1500, 1000);
			EmbeddedCanvas canXepgEphibin = new EmbeddedCanvas();
			frameXepgEphibin.add(canXepgEphibin);
			frameXepgEphibin.setLocationRelativeTo(null);
			frameXepgEphibin.setVisible(true);
			canXepgEphibin.divide(3, 2);
			for(int seci = 0; seci < 6; seci++){
				canXepgEphibin.cd(seci);
				canXepgEphibin.setFont("Arial");
				histerrGroups_X_epg_E_vs_phi.getItem(seci, theta_bini).setTitle("e^- Hit: " + theta_bnd[theta_bini] + "#degree < #theta <= "
																		+ theta_bnd[theta_bini+1] + " #degree, E_X_(ep->e'p'#gamma) vs. #phi");
				histerrGroups_X_epg_E_vs_phi.getItem(seci, theta_bini).setTitleX("#phi [#degree]");
				histerrGroups_X_epg_E_vs_phi.getItem(seci, theta_bini).setTitleY("Sector " + (seci+1) + " E_X : #mu +/- #sigma [GeV]");
				histerrGroups_X_epg_E_vs_phi.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_X_epg_E_vs_phi.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_X_epg_E_vs_phi.getItem(seci, theta_bini).setLineThickness(1);
				canXepgEphibin.getPad(seci).setTitleFontSize(32);
				canXepgEphibin.getPad(seci).setAxisTitleFontSize(32);
				canXepgEphibin.getPad(seci).setAxisLabelFontSize(24);
				canXepgEphibin.getPad(seci).setAxisRange((phi_rot[seci]-30), (phi_rot[seci]+30), -2.25, 2.25);
				histerrGroups_p_cor_X_epg_E_vs_phi.getItem(seci, theta_bini).setMarkerSize(5);
				histerrGroups_p_cor_X_epg_E_vs_phi.getItem(seci, theta_bini).setMarkerStyle(1);
				histerrGroups_p_cor_X_epg_E_vs_phi.getItem(seci, theta_bini).setMarkerColor(2);
				histerrGroups_p_cor_X_epg_E_vs_phi.getItem(seci, theta_bini).setLineColor(2);
				histerrGroups_p_cor_X_epg_E_vs_phi.getItem(seci, theta_bini).setLineThickness(1);
				canXepgEphibin.draw(histerrGroups_X_epg_E_vs_phi.getItem(seci, theta_bini), "same");
				canXepgEphibin.draw(histerrGroups_p_cor_X_epg_E_vs_phi.getItem(seci, theta_bini), "same");
			}
		}
	}
}
