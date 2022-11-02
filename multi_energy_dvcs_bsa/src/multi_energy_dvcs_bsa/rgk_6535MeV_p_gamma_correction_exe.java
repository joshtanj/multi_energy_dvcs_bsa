package multi_energy_dvcs_bsa;

import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;

import org.freehep.math.minuit.FCNBase;
import org.freehep.math.minuit.FunctionMinimum;
import org.freehep.math.minuit.MnMigrad;
import org.freehep.math.minuit.MnStrategy;
import org.freehep.math.minuit.MnUserParameters;
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

public class rgk_6535MeV_p_gamma_correction_exe {

	static HipoDataSource readerdvcs = new HipoDataSource();
	
	static H2F hdvcs_theta_vs_p_gamma = new H2F("dvcs_theta_vs_p_gamma", "dvcs_theta_vs_p_gamma", 200, 0, 8, 200, 0, 45);
	static H2F hdvcs_theta_vs_phi_gamma = new H2F("dvcs_theta_vs_phi_e", "dvcs_theta_vs_phi_e", 200, -180, 180, 200, 0, 45);
	static H1F htheta_cone_gamma = new H1F("theta_cone_gamma", "theta_cone_gamma", 100, 0, 7);
	static H1F hX_eg_m = new H1F("X_eg_m", "X_eg_m", 100, 0, 2);
	static H1F hX_epg_pt = new H1F("X_epg_pt", "X_epg_pt", 100, 0, 0.8);
	static H1F hX_epg_E = new H1F("X_epg_E", "X_epg_E", 100, -2.25, 2.25);
	static H1F hpcor_theta_cone_gamma = new H1F("pcortheta_cone_gamma", "pcortheta_cone_gamma", 100, 0, 7);
	static H1F hpcor_X_eg_m = new H1F("pcorX_eg_m", "pcorX_eg_m", 100, 0, 2);
	static H1F hpcor_X_epg_pt = new H1F("pcorX_epg_pt", "pcorX_epg_pt", 100, 0, 0.8);
	static H1F hpcor_X_epg_E = new H1F("pcorX_epg_E", "pcorX_epg_E", 100, -2.25, 2.25);
	
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
			H1F hX_eg_m_sec = new H1F("hX_eg_m_sec", "hX_eg_m_sec", 100, 0, 2);
			histGroups_X_eg_m_sec.add(hX_eg_m_sec, seci);
			H1F hp_cor_X_eg_m_sec = new H1F("hp_cor_X_eg_m_sec", "hp_cor_X_eg_m_sec", 100, 0, 2);
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
					H1F hX_eg_m_vs_phi_bin = new H1F("hX_eg_m_vs_phi_bin", "hX_eg_m_vs_phi_bin", 100, 0, 2);
					histGroups_X_eg_m_theta_phi_bin.add(hX_eg_m_vs_phi_bin, seci, thetabin, phibin);
					H1F hp_cor_X_eg_m_theta_phi_bin = new H1F("hp_cor_X_eg_m_theta_phi_bin", "hp_cor_X_eg_m_theta_phi_bin", 100, 0, 2);
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

	static void processEventdvcs(DataEvent event, int corEventCounter, IndexedList<double[]> thetaeCorFac, IndexedList<double[]> thetagammaCorFac) {
	
	if(event.hasBank("REC::Particle") && event.hasBank("REC::Track") && event.hasBank("REC::Traj"))
	{
		HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
		int ecount = 0;
		int e_pindex = -1;
		LorentzVector p_e = new LorentzVector(0, 0, 0, 0);
		LorentzVector cor_p_e = new LorentzVector(0, 0, 0, 0);
		Vector3D v_e = new Vector3D (0, 0, 0);
		float M_e = (float) 0.000511;
		double E = 6.535;
		float theta_deg_e = 0;
		float phi_deg_e = -200;
		float phi_rot_deg_e = -200;
		double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
		int protoncount = 0;
		int proton_pindex = -1;
		LorentzVector p_proton = new LorentzVector(0, 0, 0, 0);
		Vector3D v_proton = new Vector3D (0, 0, 0);
		float M_p = (float) 0.938272;
		float cdchi2_proton = 100;
		short cdNDF_proton = 100;
		int gammacount = 0;
		int gamma_pindex = -1;		
		LorentzVector p_gamma = new LorentzVector(0, 0, 0, 0);
		LorentzVector p_cor_gamma = new LorentzVector(0, 0, 0, 0);
		byte gamma_sector = 0;
		float theta_deg_gamma = 0;
		float phi_deg_gamma = -200;
		float phi_rot_deg_gamma = -200;
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
		Vector3D ecal_hit_gamma = new Vector3D (0, 0, 0);
		Vector3D ecal_hit_rot_gamma = new Vector3D (0, 0, 0);
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
						float z = reccal.getFloat("z", l);
						float lu = reccal.getFloat("lu", l);
						float lv = reccal.getFloat("lv", l);
						float lw = reccal.getFloat("lw", l);
						byte sector = reccal.getByte("sector", l);
						gamma_sector = sector;
						if(pindex == gamma_pindex && detector == 7 && layer == 1)
						{
							ecal_hit_gamma = new Vector3D (x, y, z);
							theta_deg_gamma = (float) (57.3*ecal_hit_gamma.theta());
							phi_deg_gamma = (float) (57.3*ecal_hit_gamma.phi());
							if(phi_deg_gamma <= -150) phi_deg_gamma = phi_deg_gamma+360;
							ecal_hit_rot_gamma = new Vector3D (x*Math.cos(phi_rot[sector-1]/57.3)+y*Math.sin(phi_rot[sector-1]/57.3),
															y*Math.cos(phi_rot[sector-1]/57.3)-x*Math.sin(phi_rot[sector-1]/57.3), z);
							phi_rot_deg_gamma = (float) (57.3*ecal_hit_rot_gamma.phi());
							double[] b0 = thetagammaCorFac.getItem(gamma_sector, 0);
							double[] b1 = thetagammaCorFac.getItem(gamma_sector, 1);
							
							double a0 = (b0[4]*theta_deg_gamma*theta_deg_gamma*theta_deg_gamma*theta_deg_gamma)+(b0[3]*theta_deg_gamma*theta_deg_gamma*theta_deg_gamma)
											+(b0[2]*theta_deg_gamma*theta_deg_gamma)+(b0[1]*theta_deg_gamma)+b0[0];
							double a1 = (b1[4]*theta_deg_gamma*theta_deg_gamma*theta_deg_gamma*theta_deg_gamma)+(b1[3]*theta_deg_gamma*theta_deg_gamma*theta_deg_gamma)
											+(b1[2]*theta_deg_gamma*theta_deg_gamma)+(b1[1]*theta_deg_gamma)+b1[0];
							
							ecal_hit_rot_gamma = new Vector3D (x*Math.cos(phi_rot[gamma_sector-1]/57.3)+y*Math.sin(phi_rot[gamma_sector-1]/57.3),
															y*Math.cos(phi_rot[gamma_sector-1]/57.3)-x*Math.sin(phi_rot[gamma_sector-1]/57.3), z);
							phi_rot_deg_gamma = (float) (57.3*dc_hit_rot_e.phi());
							pcf = (a1*phi_rot_deg_gamma)+a0;
							p_cor_gamma = new LorentzVector(pcf*p_gamma.px(), pcf*p_gamma.py(), pcf*p_gamma.pz(),
															Math.sqrt((pcf*p_gamma.p()*pcf*p_gamma.p())));
						}
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
				
				double[] b0 = thetaeCorFac.getItem(e_sector, 0);
				double[] b1 = thetaeCorFac.getItem(e_sector, 1);
				
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
				p_e = cor_p_e;
			}
		}
		
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
		p_cor_p_X_ep.sub(p_e);
		p_cor_p_X_ep.sub(p_proton);
		LorentzVector p_cor_p_X_eg = new LorentzVector(0, 0, 0, 0);
		p_cor_p_X_eg.add(p_B);
		p_cor_p_X_eg.add(p_T);
		p_cor_p_X_eg.sub(p_e);
		p_cor_p_X_eg.sub(p_cor_gamma);
		LorentzVector p_cor_p_X_epg = new LorentzVector(0, 0, 0, 0);
		p_cor_p_X_epg.add(p_B);
		p_cor_p_X_epg.add(p_T);
		p_cor_p_X_epg.sub(p_e);
		p_cor_p_X_epg.sub(p_proton);
		p_cor_p_X_epg.sub(p_cor_gamma);
		
		double x_B = Q2/(2*p_T.e()*q.e());
		LorentzVector t_cal = new LorentzVector(p_X_eg.px(), p_X_eg.py(), p_X_eg.pz(), p_X_eg.e());
		t_cal.sub(p_T);
		LorentzVector p_cor_t_cal = new LorentzVector(p_cor_p_X_eg.px(), p_cor_p_X_eg.py(), p_cor_p_X_eg.pz(), p_cor_p_X_eg.e());
		p_cor_t_cal.sub(p_T);
		
		if(ecount == 1 && protoncount == 1 && gammacount == 1
				&& v_e.z() > -12 && v_e.z() < 7 && Math.abs(v_e.z()-v_proton.z()) < 2.5 + (2.5/p_proton.p())
				&& p_gamma.vect().theta(p_e.vect()) > 5
				&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30))
				&& 57.3*p_proton.theta() < 75
				&& fid_ecal_gamma == true
				&& theta_deg_e > 7 && theta_deg_e <= 30
				&& Q2 > 1 && W > 2)
		{
			if(-t_cal.mass2() > 0.938272*0.938272*x_B*x_B/(1-x_B+x_B*0.938272*0.938272/Q2*Q2))
			{
				hdvcs_theta_vs_p_gamma.fill(p_gamma.p(), 57.3*p_gamma.theta());
				hdvcs_theta_vs_phi_gamma.fill(57.3*p_gamma.phi(), 57.3*p_gamma.theta());
				htheta_cone_gamma.fill(p_gamma.vect().theta(p_X_ep.vect()));
				hX_eg_m.fill(p_X_eg.mass());
				hX_epg_pt.fill(Math.sqrt(p_X_epg.px()*p_X_epg.px()+p_X_epg.py()*p_X_epg.py()));
				hX_epg_E.fill(p_X_epg.e());			
				histGroups_X_eg_m_sec.getItem(gamma_sector-1).fill(p_X_eg.mass());
			}
			
			if(-p_cor_t_cal.mass2() > 0.938272*0.938272*x_B*x_B/(1-x_B+x_B*0.938272*0.938272/Q2*Q2))
			{
					hpcor_theta_cone_gamma.fill(p_cor_gamma.vect().theta(p_cor_p_X_ep.vect()));
					hpcor_X_eg_m.fill(p_cor_p_X_eg.mass());
					hpcor_X_epg_pt.fill(Math.sqrt(p_cor_p_X_epg.px()*p_cor_p_X_epg.px()+p_cor_p_X_epg.py()*p_cor_p_X_epg.py()));
					hpcor_X_epg_E.fill(p_cor_p_X_epg.e());
					histGroups_p_cor_X_eg_m_sec.getItem(gamma_sector-1).fill(p_cor_p_X_eg.mass());
			}
			
			double[] theta_bnd  = new double[] {5, 7, 9, 11, 14, 17, 20, 25, 30, 35};
			IndexedList<double[]> phi_rot_bnd = new IndexedList<double[]>(1);
			phi_rot_bnd.add(new double[] {-30, -4, 0, 4, 30}, 0);
			phi_rot_bnd.add(new double[] {-30, -7, 0, 7, 30}, 1);
			phi_rot_bnd.add(new double[] {-30, -8, 0, 8, 30}, 2);
			phi_rot_bnd.add(new double[] {-30, -10, 0, 10, 30}, 3);
			phi_rot_bnd.add(new double[] {-30, -11, 0, 11, 30}, 4);
			phi_rot_bnd.add(new double[] {-30, -12, 0, 12, 30}, 5);
			phi_rot_bnd.add(new double[] {-30, -13, 0, 13, 30}, 6);
			phi_rot_bnd.add(new double[] {-30, -14, 0, 14, 30}, 7);
			phi_rot_bnd.add(new double[] {-30, -17, 0, 17, 30}, 8);
			for(int thbin = 0; thbin < 9; thbin++)
			{
				if(theta_deg_e > theta_bnd[thbin] && theta_deg_e <= theta_bnd[thbin+1])
				{
					for(int phbin = 0; phbin < 4; phbin++)
					{
						if(phi_rot_deg_e > phi_rot_bnd.getItem(thbin)[phbin] && phi_rot_deg_e <= phi_rot_bnd.getItem(thbin)[phbin+1])
						{
							if(-t_cal.mass2() > 0.938272*0.938272*x_B*x_B/(1-x_B+x_B*0.938272*0.938272/Q2*Q2))
							{
								histGroups_theta_cone_gamma_theta_phi_bin.getItem(gamma_sector-1, thbin, phbin)
																			.fill(p_gamma.vect().theta(p_X_ep.vect()));	
								histGroups_X_eg_m_theta_phi_bin.getItem(gamma_sector-1, thbin, phbin).fill(p_X_eg.mass());	
								histGroups_X_epg_pt_theta_phi_bin.getItem(gamma_sector-1, thbin, phbin)
																	.fill(Math.sqrt(p_X_epg.px()*p_X_epg.px()+p_X_epg.py()*p_X_epg.py()));		
								histGroups_X_epg_E_theta_phi_bin.getItem(gamma_sector-1, thbin, phbin).fill(p_X_epg.e());
							}
							if(-p_cor_t_cal.mass2() > 0.938272*0.938272*x_B*x_B/(1-x_B+x_B*0.938272*0.938272/Q2*Q2))
							{
								histGroups_p_cor_theta_cone_gamma_theta_phi_bin.getItem(gamma_sector-1, thbin, phbin)
																				.fill(p_gamma.vect().theta(p_cor_p_X_ep.vect()));
								histGroups_p_cor_X_eg_m_theta_phi_bin.getItem(gamma_sector-1, thbin, phbin).fill(p_cor_p_X_eg.mass());
								histGroups_p_cor_X_epg_pt_theta_phi_bin.getItem(gamma_sector-1, thbin, phbin)
																		.fill(Math.sqrt(p_cor_p_X_epg.px()*p_cor_p_X_epg.px()
																						+p_cor_p_X_epg.py()*p_cor_p_X_epg.py()));
								histGroups_p_cor_X_epg_E_theta_phi_bin.getItem(gamma_sector-1, thbin, phbin).fill(p_cor_p_X_epg.e());
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
	
		histos_dvcs_ev_sec();
		histos_dvcs_ev_vs_phi_bin();
		
		IndexedList<double[]> thetaeCorFac = new IndexedList<>(2);
		double[] thetaeCorFac_sec1_b0 = new double[] {0.9938077481295297, 8.258708800790782E-4, -2.2658532616476987E-5,
														9.677115140791733E-7, -4.5642462769384737E-8};
		double[] thetaeCorFac_sec1_b1 = new double[] {0.0024178549687052706, -6.653055638986705E-4, 7.219407837197144E-5,
														-3.180462703382774E-6, 4.9352098479104826E-8};
		double[] thetaeCorFac_sec2_b0 = new double[] {1.0035661097100197, -0.003285642462705683, 5.292556662659571E-4,
														-3.0231148934779212E-5, 5.609421915775323E-7};
		double[] thetaeCorFac_sec2_b1 = new double[] {5.346964410100895E-4, -6.714148268661577E-5, 5.201475751662365E-6,
														-1.6155087929089417E-7, 1.0186220422620223E-9};
		double[] thetaeCorFac_sec3_b0 = new double[] {1.0087429705449467, -0.004408627402638188, 6.482030554858202E-4,
														-3.394142844711124E-5, 5.779315008382153E-7};
		double[] thetaeCorFac_sec3_b1 = new double[] {-4.397136504767579E-4, 1.4884197765913248E-4, -1.8433837666625867E-5,
														8.244944205916959E-7, -1.228942642526823E-8};
		double[] thetaeCorFac_sec4_b0 = new double[] {1.006956918896026, -0.003711867838701133, 5.518086095685475E-4,
														-2.8746618476848303E-5, 4.828563256249679E-7};
		double[] thetaeCorFac_sec4_b1 = new double[] {8.309626628189099E-4, -2.701812621232045E-4, 2.875524544466852E-5,
														-1.2835027495862537E-6, 2.029581677685505E-8};
		double[] thetaeCorFac_sec5_b0 = new double[] {0.9570390799595814, 0.009328675535642626, -6.714343815803992E-4,
														1.987029852607302E-5, -2.1755105045221692E-7};
		double[] thetaeCorFac_sec5_b1 = new double[] {-5.873762378695643E-4, 1.9121282318863288E-4, -1.9101091097978772E-5,
														7.82569493538631E-7, -1.1798923261435443E-8};
		double[] thetaeCorFac_sec6_b0 = new double[] {0.9636208140401962, 0.00948121762243233, -7.335233524854341E-4,
														2.4169978390757133E-5, -3.071125516271284E-7};
		double[] thetaeCorFac_sec6_b1 = new double[] {0.0017450577912736217, -5.245913063775676E-4, 5.3940149514332456E-5,
														-2.3107384729861627E-6, 3.575938540124231E-8};
		thetaeCorFac.add(thetaeCorFac_sec1_b0, 1, 0);
		thetaeCorFac.add(thetaeCorFac_sec1_b1, 1, 1);
		thetaeCorFac.add(thetaeCorFac_sec2_b0, 2, 0);
		thetaeCorFac.add(thetaeCorFac_sec2_b1, 2, 1);
		thetaeCorFac.add(thetaeCorFac_sec3_b0, 3, 0);
		thetaeCorFac.add(thetaeCorFac_sec3_b1, 3, 1);
		thetaeCorFac.add(thetaeCorFac_sec4_b0, 4, 0);
		thetaeCorFac.add(thetaeCorFac_sec4_b1, 4, 1);
		thetaeCorFac.add(thetaeCorFac_sec5_b0, 5, 0);
		thetaeCorFac.add(thetaeCorFac_sec5_b1, 5, 1);
		thetaeCorFac.add(thetaeCorFac_sec6_b0, 6, 0);
		thetaeCorFac.add(thetaeCorFac_sec6_b1, 6, 1);
		
		double[] theta_bnd  = new double[] {5, 7, 9, 11, 14, 17, 20, 25, 30, 35};
		double[] theta_ba = new double[]{6.28, 7.99, 9.97, 12.48, 15.47, 18.47, 22.41, 27.46, 31.45};
		double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
		IndexedList<double[]> phi_rot_bnd = new IndexedList<double[]>(1);
		phi_rot_bnd.add(new double[] {-30, -4, 0, 4, 30}, 0);
		phi_rot_bnd.add(new double[] {-30, -7, 0, 7, 30}, 1);
		phi_rot_bnd.add(new double[] {-30, -8, 0, 8, 30}, 2);
		phi_rot_bnd.add(new double[] {-30, -10, 0, 10, 30}, 3);
		phi_rot_bnd.add(new double[] {-30, -11, 0, 11, 30}, 4);
		phi_rot_bnd.add(new double[] {-30, -12, 0, 12, 30}, 5);
		phi_rot_bnd.add(new double[] {-30, -13, 0, 13, 30}, 6);
		phi_rot_bnd.add(new double[] {-30, -14, 0, 14, 30}, 7);
		phi_rot_bnd.add(new double[] {-30, -17, 0, 17, 30}, 8);
		IndexedList<double[]> phi_rot_ba = new IndexedList<double[]>(1);
		phi_rot_ba.add(new double[] {-7.36, -2.24, 2.24, 7.36}, 0);
		phi_rot_ba.add(new double[] {-10.46, -3.45, 3.45, 10.46}, 1);
		phi_rot_ba.add(new double[] {-12.74, -4.18, 4.18, 12.74}, 2);
		phi_rot_ba.add(new double[] {-14.62, -4.97, 4.97, 14.62}, 3);
		phi_rot_ba.add(new double[] {-16.17, -5.67, 5.67, 16.17}, 4);
		phi_rot_ba.add(new double[] {-17.44, -6.30, 6.30, 17.44}, 5);
		phi_rot_ba.add(new double[] {-18.48, -6.90, 6.90, 18.48}, 6);
		phi_rot_ba.add(new double[] {-19.23, -7.20, 7.20, 19.23}, 7);
		phi_rot_ba.add(new double[] {-21.35, -9.42, 9.42, 21.35}, 8);
		
		IndexedList<double[]> thetagammaCorFac = new IndexedList<>(2);
		double[] thetagammaCorFac_sec1_b0 = new double[] {1.1243869657347507, -0.017306168428098012, 0.0015898472292506453, 
				-3.941806978097773E-5, 3.266397723426338E-7};
double[] thetagammaCorFac_sec1_b1 = new double[] {-0.004976450843098421, 0.00131383666530172, -1.2249208744308194E-4, 
				4.521856751725047E-6, -5.711534486885231E-8};
double[] thetagammaCorFac_sec2_b0 = new double[] {1.3176609850916636, -0.08480658598284278, 0.008475310072293078, 
				-3.0796169213389505E-4, 3.913584565103736E-6};
double[] thetagammaCorFac_sec2_b1 = new double[] {0.0065829400132483805, -0.0018361428092863556, 1.6372584633989049E-4, 
				-6.139343249964539E-6, 8.251353447607186E-8};
double[] thetagammaCorFac_sec3_b0 = new double[] {1.4160632427178463, -0.10907935076383407, 0.01062757653019337, 
				-3.8956339738527876E-4, 5.005131959548998E-6};
double[] thetagammaCorFac_sec3_b1 = new double[] {-0.0016375401054715744, 1.4838699180376688E-4, -4.809863943880121E-6, 
				-1.3191921487381088E-7, 6.481916186027909E-9};
double[] thetagammaCorFac_sec4_b0 = new double[] {1.1422778694371523, -0.033330975681096914, 0.00393741516411168, 
				-1.4680852848429277E-4, 1.8533360277899754E-6};
double[] thetagammaCorFac_sec4_b1 = new double[] {-0.0028583797769678866, 8.10171716027917E-4, -8.910617277522397E-5, 
				3.6165878224665257E-6, -4.843896846702168E-8};
double[] thetagammaCorFac_sec5_b0 = new double[] {1.2867713580270634, -0.07227527249118218, 0.007307399954943301, 
				-2.62975386884705E-4, 3.278091096321523E-6};
double[] thetagammaCorFac_sec5_b1 = new double[] {0.0034551337942566145, -0.0011080141834231391, 9.938654253057621E-5, 
				-3.7222476021270467E-6, 4.986981387012025E-8};
double[] thetagammaCorFac_sec6_b0 = new double[] {1.3707719748308587, -0.09194115239209903, 0.008895676243258127, 
				-3.184343551517498E-4, 3.982133189812425E-6};
double[] thetagammaCorFac_sec6_b1 = new double[] {0.005179116815909332, -0.0016018753298132156, 1.5713477874125626E-4, 
				-6.2713876479417714E-6, 8.771089167411175E-8};

		thetagammaCorFac.add(thetagammaCorFac_sec1_b0, 1, 0);
		thetagammaCorFac.add(thetagammaCorFac_sec1_b0, 1, 0);
		thetagammaCorFac.add(thetagammaCorFac_sec1_b1, 1, 1);
		thetagammaCorFac.add(thetagammaCorFac_sec2_b0, 2, 0);
		thetagammaCorFac.add(thetagammaCorFac_sec2_b1, 2, 1);
		thetagammaCorFac.add(thetagammaCorFac_sec3_b0, 3, 0);
		thetagammaCorFac.add(thetagammaCorFac_sec3_b1, 3, 1);
		thetagammaCorFac.add(thetagammaCorFac_sec4_b0, 4, 0);
		thetagammaCorFac.add(thetagammaCorFac_sec4_b1, 4, 1);
		thetagammaCorFac.add(thetagammaCorFac_sec5_b0, 5, 0);
		thetagammaCorFac.add(thetagammaCorFac_sec5_b1, 5, 1);
		thetagammaCorFac.add(thetagammaCorFac_sec6_b0, 6, 0);
		thetagammaCorFac.add(thetagammaCorFac_sec6_b1, 6, 1);
		
	// DVCS	
		
		readerdvcs.open("C:/Users/joshtanj/Documents/download/skim_6bank_epg_merged_6535MeV_skim16.hipo");
		
		int dvcsEventCounter = 0;
		while(readerdvcs.hasEvent())// && dvcsEventCounter < 10000000)
		{
			dvcsEventCounter++;
			processEventdvcs(readerdvcs.getNextEvent(), dvcsEventCounter, thetaeCorFac, thetagammaCorFac);
			if(dvcsEventCounter%500000 == 0) System.out.println("DVCS Event: " + dvcsEventCounter);
		}

		readerdvcs.close();
		
		F1D ftc;
		F1D fpm;
		F1D fXpt;
		F1D fXe;
	
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
												.addPoint((phi_rot_ba.getItem(theta_bini)[phi_bini]+phi_rot[seci]),
															meanthetaconegamma.getItem(seci, theta_bini,  phi_bini),
															0, sigmathetaconegamma.getItem(seci, theta_bini, phi_bini));
					if(meanpcorthetaconegamma.hasItem(seci, theta_bini,  phi_bini) && sigmapcorthetaconegamma.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_p_cor_theta_cone_gamma_vs_phi.getItem(seci, theta_bini)
												.addPoint((phi_rot_ba.getItem(theta_bini)[phi_bini]+phi_rot[seci]),
															meanpcorthetaconegamma.getItem(seci, theta_bini,  phi_bini),
															0, sigmapcorthetaconegamma.getItem(seci, theta_bini, phi_bini));
					
					if(meanXegm.hasItem(seci, theta_bini,  phi_bini) && sigmaXegm.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_X_eg_m_vs_phi.getItem(seci, theta_bini)
												.addPoint((phi_rot_ba.getItem(theta_bini)[phi_bini]+phi_rot[seci]),
															meanXegm.getItem(seci, theta_bini,  phi_bini),
															0, sigmaXegm.getItem(seci, theta_bini, phi_bini));
					if(meanpcorXegm.hasItem(seci, theta_bini,  phi_bini) && sigmapcorXegm.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_p_cor_X_eg_m_vs_phi.getItem(seci, theta_bini)
												.addPoint((phi_rot_ba.getItem(theta_bini)[phi_bini]+phi_rot[seci]),
															meanpcorXegm.getItem(seci, theta_bini,  phi_bini),
															0, sigmapcorXegm.getItem(seci, theta_bini, phi_bini));
					
					if(meanXepgpt.hasItem(seci, theta_bini,  phi_bini) && sigmaXepgpt.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_X_epg_pt_vs_phi.getItem(seci, theta_bini)
												.addPoint((phi_rot_ba.getItem(theta_bini)[phi_bini]+phi_rot[seci]),
															meanXepgpt.getItem(seci, theta_bini,  phi_bini),
															0, sigmaXepgpt.getItem(seci, theta_bini, phi_bini));
					if(meanpcorXepgpt.hasItem(seci, theta_bini,  phi_bini) && sigmapcorXepgpt.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_p_cor_X_epg_pt_vs_phi.getItem(seci, theta_bini)
												.addPoint((phi_rot_ba.getItem(theta_bini)[phi_bini]+phi_rot[seci]),
															meanpcorXepgpt.getItem(seci, theta_bini,  phi_bini),
															0, sigmapcorXepgpt.getItem(seci, theta_bini, phi_bini));
					
					if(meanXepgE.hasItem(seci, theta_bini,  phi_bini) && sigmaXepgE.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_X_epg_E_vs_phi.getItem(seci, theta_bini)
												.addPoint((phi_rot_ba.getItem(theta_bini)[phi_bini]+phi_rot[seci]),
															meanXepgE.getItem(seci, theta_bini,  phi_bini),
															0, sigmaXepgE.getItem(seci, theta_bini, phi_bini));
					if(meanpcorXepgE.hasItem(seci, theta_bini,  phi_bini) && sigmapcorXepgE.hasItem(seci, theta_bini,  phi_bini))
						histerrGroups_p_cor_X_epg_E_vs_phi.getItem(seci, theta_bini)
												.addPoint((phi_rot_ba.getItem(theta_bini)[phi_bini]+phi_rot[seci]),
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
