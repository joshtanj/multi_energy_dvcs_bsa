package multi_energy_dvcs_bsa;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.swing.JFrame;

import org.jlab.clas.physics.LorentzVector;
import org.jlab.geom.prim.Vector3D;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.utils.groups.IndexedList;
import org.slf4j.Marker;
import org.jlab.io.evio.EvioDataBank;
import org.jlab.io.hipo.HipoDataBank;

public class rgk_6535MeV_dvcs_analysis_bsa_sub_pi0_sim_var_tbin {
	
	static HipoDataSource reader = new HipoDataSource();
	static HipoDataSource reader_pi0 = new HipoDataSource();
	static HipoDataSource reader_pi0_sim = new HipoDataSource();
	
	static H1F hphi_poshel = new H1F("phi_poshel", "phi_poshel", 120, 0, 360);
	static H1F hphi_neghel = new H1F("phi_neghel", "phi_neghel", 120, 0, 360);
	static H1F hphi_poshel_pi0 = new H1F("phi_poshel_pi0", "phi_poshel_pi0", 120, 0, 360);
	static H1F hphi_neghel_pi0 = new H1F("phi_neghel_pi0", "phi_neghel_pi0", 120, 0, 360);
	static H1F hphi_1gamma_pi0_sim = new H1F("phi_1gamma_pi0_sim", "phi_1gamma_pi0_sim", 120, 0, 360);
	static H1F hphi_2gamma_pi0_sim = new H1F("phi_2gamma_pi0_sim", "phi_2gamma_pi0_sim", 120, 0, 360);
	
	static IndexedList<H1F> histGroups_phi_poshel_tbin = new IndexedList<H1F>(2);
	static IndexedList<H1F> histGroups_phi_poshel_tbin_pi0 = new IndexedList<H1F>(2);
	public static void phi_poshel_tbin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			for(int itbin = 0; itbin < 5; itbin++)
			{
				H1F phi_poshel_tbin = new H1F("phi_poshel_tbin", "phi_poshel_tbin", 120, 0, 360);
				histGroups_phi_poshel_tbin.add(phi_poshel_tbin, ibin, itbin);
				H1F phi_poshel_tbin_pi0 = new H1F("phi_poshel_tbin_pi0", "phi_poshel_tbin_pi0", 120, 0, 360);
				histGroups_phi_poshel_tbin_pi0.add(phi_poshel_tbin_pi0, ibin, itbin);
			}
		}
	}
	
	static IndexedList<H1F> histGroups_phi_neghel_tbin = new IndexedList<H1F>(2);
	static IndexedList<H1F> histGroups_phi_neghel_tbin_pi0 = new IndexedList<H1F>(2);
	public static void phi_neghel_tbin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			for(int itbin = 0; itbin < 5; itbin++)
			{
				H1F phi_neghel_tbin = new H1F("phi_neghel_tbin", "phi_neghel_tbin", 120, 0, 360);
				histGroups_phi_neghel_tbin.add(phi_neghel_tbin, ibin, itbin);
				H1F phi_neghel_tbin_pi0 = new H1F("phi_neghel_tbin_pi0", "phi_neghel_tbin_pi0", 120, 0, 360);
				histGroups_phi_neghel_tbin_pi0.add(phi_neghel_tbin_pi0, ibin, itbin);
			}
		}
	}
	
	static IndexedList<H1F> histGroups_phi_1gamma_pi0_sim_tbin = new IndexedList<H1F>(2);
	public static void phi_1gamma_pi0_sim_tbin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			for(int itbin = 0; itbin < 5; itbin++)
			{
				H1F phi_1gamma_pi0_sim_tbin = new H1F("phi_1gamma_pi0_sim_tbin", "phi_1gamma_pi0_sim_tbin", 120, 0, 360);
				histGroups_phi_1gamma_pi0_sim_tbin.add(phi_1gamma_pi0_sim_tbin, ibin, itbin);
			}
		}
	}
	
	static IndexedList<H1F> histGroups_phi_2gamma_pi0_sim_tbin = new IndexedList<H1F>(2);
	public static void phi_2gamma_pi0_sim_tbin_histos() {
		for(int ibin = 0; ibin  < 8; ibin++) {
			for(int itbin = 0; itbin < 5; itbin++)
			{
				H1F phi_2gamma_pi0_sim_tbin = new H1F("phi_2gamma_pi0_sim_tbin", "phi_2gamma_pi0_sim_tbin", 120, 0, 360);
				histGroups_phi_2gamma_pi0_sim_tbin.add(phi_2gamma_pi0_sim_tbin, ibin, itbin);
			}
		}
	}
	
	static double rc = Math.PI/180;
	
	static void processEvent(DataEvent event, int eventCounter)
	{
		if(event.hasBank("REC::Particle") && event.hasBank("REC::Track"))
		{
			HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
			int ecount = 0;
			int protoncount = 0;
			int gammacount = 0;
			int e_pindex = -1;
			boolean fid_ecal_e = false;
			int proton_pindex = -1;
			int gamma1_pindex = -1;
			int gamma2_pindex = -1;
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0.000511);
			LorentzVector p_proton = new LorentzVector(0, 0, 0, 0.938272);
			LorentzVector p_gamma1 = new LorentzVector(0, 0, 0, 0);
			LorentzVector p_gamma2 = new LorentzVector(0, 0, 0, 0);
			double E = 6.535;
			double E_gamma_cut = 0.3;
			byte det_gamma1 = 0;
			byte det_gamma2 = 0;
			LorentzVector v_e = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_proton = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_gamma1 = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_gamma2 = new LorentzVector(0, 0, 0, 0);
			float cdchi2_proton = 100;
			short cdNDF_proton = 100;
			double cX = 0;
			double cY = 0;
			boolean fid_ecal_gamma1 = false;
			boolean fid_ecal_gamma2 = false;
			boolean fid_ftcal_gamma1 = false;
			boolean fid_ftcal_gamma2 = false;
			double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
			double[] cX_split  = new double[] {87, 82, 85, 77, 78, 82};
			double[] t_left  = new double[] {58.7356, 62.8204, 62.2296, 53.7756, 58.2888, 54.5822};
			double[] t_right  = new double[] {58.7477, 51.2589, 59.2357, 56.2415, 60.8219, 49.8914};
			double[] s_left  = new double[] {0.582053, 0.544976, 0.549788, 0.56899, 0.56414, 0.57343};
			double[] s_right  = new double[] {-0.591876, -0.562926, -0.562246, -0.563726, -0.568902, -0.550729};
			double[] r_left  = new double[] {64.9348, 64.7541, 67.832, 55.9324, 55.9225, 60.0997};
			double[] r_right  = new double[] {65.424, 54.6992, 63.6628, 57.8931, 56.5367, 56.4641};
			double[] q_left  = new double[] {0.745578, 0.606081, 0.729202, 0.627239, 0.503674, 0.717899};
			double[] q_right  = new double[] {-0.775022, -0.633863, -0.678901, -0.612458, -0.455319, -0.692481};
			for(int i = 0; i < rec.rows(); i++)
			{
				int pid = rec.getInt("pid", i);
				float px = rec.getFloat("px", i);
				float py = rec.getFloat("py", i);
				float pz = rec.getFloat("pz", i);
				float vx = rec.getFloat("vx", i);
				float vy = rec.getFloat("vy", i);
				float vz = rec.getFloat("vz", i);
				byte charge = rec.getByte("charge", i);
				float beta = rec.getFloat("beta", i);
				float p = (float) Math.sqrt((px*px)+(py*py)+(pz*pz));
				if(pid == 11)
				{
					ecount++;
					if(ecount == 1)
					{
						p_e = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.000511*0.000511)));
						v_e = new LorentzVector(vx, vy, vz, 0);
						e_pindex = i;
					}
				}
				if(pid == 2212)
				{
					protoncount++;
					if(protoncount == 1)
					{
						p_proton = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.938272*0.938272)));
						v_proton = new LorentzVector(vx, vy, vz, 0);
						proton_pindex = i;
					}
				}
				if(pid == 22)
				{
					gammacount++;
					if(gammacount == 1)
					{	
						p_gamma1 = new LorentzVector(px, py, pz, p);
						v_gamma1 = new LorentzVector(vx, vy, vz, 0);
						gamma1_pindex = i;
					}
					if(gammacount == 2)
					{
						p_gamma2 = new LorentzVector(px, py, pz, p);
						v_gamma2 = new LorentzVector(vx, vy, vz, 0);
						gamma2_pindex = i;
					}
				}
			}
			//Start: cal fiducial cut
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
					if(pindex == e_pindex && detector == 7 && layer == 1)
					{
						cX = x*Math.cos(rc*phi_rot[sector-1])+y*Math.sin(rc*phi_rot[sector-1]);
						cY = y*Math.cos(rc*phi_rot[sector-1])-x*Math.sin(rc*phi_rot[sector-1]);
						if(cX > cX_split[sector-1] && cY < s_left[sector-1]*(cX-t_left[sector-1])
							&& cY > s_right[sector-1]*(cX-t_right[sector-1]))
							fid_ecal_e = true;
						else if(cX < cX_split[sector-1] && cY < q_left[sector-1]*(cX-r_left[sector-1])
									&& cY > q_right[sector-1]*(cX-r_right[sector-1]))
							fid_ecal_e = true;
					}
					
					if(pindex == gamma1_pindex && detector == 7 && layer == 1)
					{
						cX = x*Math.cos(rc*phi_rot[sector-1])+y*Math.sin(rc*phi_rot[sector-1]);
						cY = y*Math.cos(rc*phi_rot[sector-1])-x*Math.sin(rc*phi_rot[sector-1]);
						if(cX > cX_split[sector-1] && cY < s_left[sector-1]*(cX-t_left[sector-1])
							&& cY > s_right[sector-1]*(cX-t_right[sector-1]))
							fid_ecal_gamma1 = true;
						else if(cX < cX_split[sector-1] && cY < q_left[sector-1]*(cX-r_left[sector-1])
									&& cY > q_right[sector-1]*(cX-r_right[sector-1]))
							fid_ecal_gamma1 = true;
						det_gamma1 = 7;
					}
					
					if(pindex == gamma1_pindex && detector == 7 && sector == 1)
					{
						if(layer == 1)
						{
							if((lw > 74.5 && lw < 79.5) || (lw > 85.5 && lw < 90.5) || (lw > 213 && lw < 218)
									|| (lw > 224.5 && lw < 229.5) || (lw > 410)) fid_ecal_gamma1 = false;
						}
						else if(layer == 4)
						{
							if((lv > 70 && lv < 93)) fid_ecal_gamma1 = false;
						}
						else if(layer == 7)
						{
							if((lu > 410.5)) fid_ecal_gamma1 = false;
							if((lv < 24) || (lv > 36.5 && lv < 60) || (lv > 411)) fid_ecal_gamma1 = false;
							if((lw < 21.5)) fid_ecal_gamma1 = false;
						}
					}
					else if(pindex == gamma1_pindex && detector == 7 && sector == 2)
					{
						if(layer == 1)
						{
							if((lv > 102 && lv < 113)) fid_ecal_gamma1 = false;
						}
						else if(layer == 4)
						{
							if((lu > 396)) fid_ecal_gamma1 = false;
							if((lw > 363)) fid_ecal_gamma1 = false;
						}
						else if(layer == 7)
						{
							if((lu < 12)) fid_ecal_gamma1 = false;
							if((lw < 10.5) || (lw > 376)) fid_ecal_gamma1 = false;
						}
					}
					else if(pindex == gamma1_pindex && detector == 7 && sector == 3)
					{
						if(layer == 1)
						{
							if((lv > 354.5 && lv < 376.5)) fid_ecal_gamma1 = false;
						}
						else if(layer == 4)
						{
							if((lu < 23)) fid_ecal_gamma1 = false;
							if((lw < 10) || (lw > 363)) fid_ecal_gamma1 = false;
						}
						else if(layer == 7)
						{
							if((lw > 387)) fid_ecal_gamma1 = false;
						}
					}
					else if(pindex == gamma1_pindex && detector == 7 && sector == 4)
					{
						if(layer == 1)
						{
							if((lv < 13) || (lv > 230 && lv < 240.5)) fid_ecal_gamma1 = false;
							if((lw > 410)) fid_ecal_gamma1 = false;
						}
						else if(layer == 4)
						{
							if((lu < 20.5)) fid_ecal_gamma1 = false;
						}
						else if(layer == 7)
						{
							if((lw < 32.5)) fid_ecal_gamma1 = false;
						}
					}
					else if(pindex == gamma1_pindex && detector == 7 && sector == 5)
					{
						if(layer == 4)
						{
							if((lv < 23)) fid_ecal_gamma1 = false;
							if((lw < 10)) fid_ecal_gamma1 = false;
						}
						else if(layer == 7)
						{
							if((lu > 193.5 && lu < 217)) fid_ecal_gamma1 = false;
							if((lv < 24)) fid_ecal_gamma1 = false;
						}
					}
					else if(pindex == gamma1_pindex && detector == 7 && sector == 6)
					{
						if(layer == 1)
						{
							if((lw > 174.5 && lw < 179.5) || (lw > 185.5 && lw < 190.5)) fid_ecal_gamma1 = false;
						}
						else if(layer == 4)
						{
							if((lv < 11.5)) fid_ecal_gamma1 = false;
							if((lu < 20.5)) fid_ecal_gamma1 = false;
						}
						else if(layer == 7)
						{
							if((lv < 12) || (lv > 423)) fid_ecal_gamma1 = false;
							if((lw < 32.5)) fid_ecal_gamma1 = false;
						}
					}

					if(pindex == gamma2_pindex && detector == 7 && layer == 1)
					{
						cX = x*Math.cos(rc*phi_rot[sector-1])+y*Math.sin(rc*phi_rot[sector-1]);
						cY = y*Math.cos(rc*phi_rot[sector-1])-x*Math.sin(rc*phi_rot[sector-1]);
						if(cX > cX_split[sector-1] && cY < s_left[sector-1]*(cX-t_left[sector-1])
							&& cY > s_right[sector-1]*(cX-t_right[sector-1]))
							fid_ecal_gamma2 = true;
						else if(cX < cX_split[sector-1] && cY < q_left[sector-1]*(cX-r_left[sector-1])
									&& cY > q_right[sector-1]*(cX-r_right[sector-1]))
							fid_ecal_gamma2 = true;
							det_gamma2 = 7;
					}
					
					if(pindex == gamma2_pindex && detector == 7 && sector == 1)
					{
						if(layer == 1)
						{
							if((lw > 74.5 && lw < 79.5) || (lw > 85.5 && lw < 90.5) || (lw > 213 && lw < 218)
									|| (lw > 224.5 && lw < 229.5) || (lw > 410)) fid_ecal_gamma2 = false;
						}
						else if(layer == 4)
						{
							if((lv > 70 && lv < 93)) fid_ecal_gamma2 = false;
						}
						else if(layer == 7)
						{
							if((lu > 410.5)) fid_ecal_gamma2 = false;
							if((lv < 24) || (lv > 36.5 && lv < 60) || (lv > 411)) fid_ecal_gamma2 = false;
							if((lw < 21.5)) fid_ecal_gamma2 = false;
						}
					}
					else if(pindex == gamma2_pindex && detector == 7 && sector == 2)
					{
						if(layer == 1)
						{
							if((lv > 102 && lv < 113)) fid_ecal_gamma2 = false;
						}
						else if(layer == 4)
						{
							if((lu > 396)) fid_ecal_gamma2 = false;
							if((lw > 363)) fid_ecal_gamma2 = false;
						}
						else if(layer == 7)
						{
							if((lu < 12)) fid_ecal_gamma2 = false;
							if((lw < 10.5) || (lw > 376)) fid_ecal_gamma2 = false;
						}
					}
					else if(pindex == gamma2_pindex && detector == 7 && sector == 3)
					{
						if(layer == 1)
						{
							if((lv > 354.5 && lv < 376.5)) fid_ecal_gamma2 = false;
						}
						else if(layer == 4)
						{
							if((lu < 23)) fid_ecal_gamma2 = false;
							if((lw < 10) || (lw > 363)) fid_ecal_gamma2 = false;
						}
						else if(layer == 7)
						{
							if((lw > 387)) fid_ecal_gamma2 = false;
						}
					}
					else if(pindex == gamma2_pindex && detector == 7 && sector == 4)
					{
						if(layer == 1)
						{
							if((lv < 13) || (lv > 230 && lv < 240.5)) fid_ecal_gamma2 = false;
							if((lw > 410)) fid_ecal_gamma2 = false;
						}
						else if(layer == 4)
						{
							if((lu < 20.5)) fid_ecal_gamma2 = false;
						}
						else if(layer == 7)
						{
							if((lw < 32.5)) fid_ecal_gamma2 = false;
						}
					}
					else if(pindex == gamma2_pindex && detector == 7 && sector == 5)
					{
						if(layer == 4)
						{
							if((lv < 23)) fid_ecal_gamma2 = false;
							if((lw < 10)) fid_ecal_gamma2 = false;
						}
						else if(layer == 7)
						{
							if((lu > 193.5 && lu < 217)) fid_ecal_gamma2 = false;
							if((lv < 24)) fid_ecal_gamma2 = false;
						}
					}
					else if(pindex == gamma2_pindex && detector == 7 && sector == 6)
					{
						if(layer == 1)
						{
							if((lw > 174.5 && lw < 179.5) || (lw > 185.5 && lw < 190.5)) fid_ecal_gamma2 = false;
						}
						else if(layer == 4)
						{
							if((lv < 11.5)) fid_ecal_gamma2 = false;
							if((lu < 20.5)) fid_ecal_gamma2 = false;
						}
						else if(layer == 7)
						{
							if((lv < 12) || (lv > 423)) fid_ecal_gamma2 = false;
							if((lw < 32.5)) fid_ecal_gamma2 = false;
						}
					}		
				}
			}
			
			if(event.hasBank("REC::ForwardTagger"))
			{
				HipoDataBank recft = (HipoDataBank) event.getBank("REC::ForwardTagger");
				for(int m = 0; m < recft.rows(); m++)
				{
					short pindex = recft.getShort("pindex", m);
					byte detector = recft.getByte("detector", m);
					byte layer = recft.getByte("layer", m);
					float x = recft.getFloat("x", m);
					float y = recft.getFloat("y", m);
					float z = recft.getFloat("z", m);
					Vector3D ftcal_det = new Vector3D (x, y, z);
					double xtrans = x + 18.36;
					double ytrans = y + 18.36;
					int cx = (int) (xtrans/1.53);
					int cy = (int) (ytrans/1.53);
					if(pindex == gamma1_pindex && detector == 10 && layer == 1)
					{
						if(ftcal_det.rho() > 8.8 && !(ftcal_det.rho() > 15.5)
							&& (x+8.5)*(x+8.5)+(y-10)*(y-10) > 1.5*1.5 
							&& (x+6)*(x+6)/(2.3*2.3)+(y+13)*(y+13)/(1.6*1.6) > 1
							&& (x-3.8)*(x-3.8)+(y+6.8)*(y+6.8) > 1.7*1.7) fid_ftcal_gamma1 = true;
						if(cx == 7 && cy == 15) fid_ftcal_gamma1 = false;
						det_gamma1 = 10;
					}
					if(pindex == gamma2_pindex && detector == 10 && layer == 1)
					{
						if(ftcal_det.rho() > 8.8 && !(ftcal_det.rho() > 15.5)
							&& (x+8.5)*(x+8.5)+(y-10)*(y-10) > 1.5*1.5 
							&& (x+6)*(x+6)/(2.3*2.3)+(y+13)*(y+13)/(1.6*1.6) > 1
							&& (x-3.8)*(x-3.8)+(y+6.8)*(y+6.8) > 1.7*1.7) fid_ftcal_gamma2 = true;
						if(cx == 7 && cy == 15) fid_ftcal_gamma2 = false;
						det_gamma2 = 10;
					}
				}
			}
			// End: cal fiducial cut
			HipoDataBank rectrac = (HipoDataBank) event.getBank("REC::Track");
			byte e_detector = 0;
			byte e_sector = 0;
			byte proton_detector = 0;
			byte gamma_sector = 0;
			byte X0_sector = 0;
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
			HipoDataBank recev = (HipoDataBank) event.getBank("REC::Event");
			byte helic = recev.getByte("helicity", 0);		
			LorentzVector q = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			LorentzVector qr = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			q.sub(p_e);
			qr.sub(p_e);
			double Q2 = -1*q.mass2();
			LorentzVector X = new LorentzVector(0, 0, 0, 0.938272);
			double x_B = Q2/(2*X.e()*q.e());
			LorentzVector t = new LorentzVector(p_proton.px(), p_proton.py(), p_proton.pz(), p_proton.e());
			t.sub(X);
			X.add(q);
			q.sub(p_gamma1);
			double W = X.mass();
			X.sub(p_proton);
			X.sub(p_gamma1);
			LorentzVector X_e = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_e.sub(p_proton);
			X_e.sub(p_gamma1);
			LorentzVector X_proton = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_proton.sub(p_e);
			X_proton.sub(p_gamma1);
			LorentzVector X_tar = new LorentzVector(0, 0, 0, 0.938272);
			LorentzVector t_cal = new LorentzVector(X_proton.px(), X_proton.py(), X_proton.pz(), X_proton.e());
			t_cal.sub(X_tar);
			LorentzVector X_gamma1 = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_gamma1.sub(p_proton);
			X_gamma1.sub(p_e);
			LorentzVector beam = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			
			double trento = ((beam.vect().cross(p_e.vect())).dot(p_proton.vect()))/Math.abs((beam.vect().cross(p_e.vect())).dot(p_proton.vect()));
			double phi = trento*57.3*Math.acos(((beam.vect().cross(p_e.vect())).dot((p_proton.vect().cross(qr.vect()))))
							/((beam.vect().cross(p_e.vect()).mag()*(p_proton.vect().cross(qr.vect())).mag())));
			if(phi < 0) phi = phi+360;
			
			int Q2xBbin = -1;
			if(x_B < 0.16) Q2xBbin = 0;
			else if(x_B >= 0.16 && x_B < 0.21 && Q2 < 1.2) Q2xBbin = 1;
			else if(x_B >= 0.16 && x_B < 0.21 && Q2 >= 1.2) Q2xBbin = 2;
			else if(x_B >= 0.21 && x_B < 0.28 && Q2 < 1.2) Q2xBbin = 3;
			else if(x_B >= 0.21 && x_B < 0.28 && Q2 >= 1.2 && Q2 < 1.6) Q2xBbin = 4;
			else if(x_B >= 0.21 && x_B < 0.28 && Q2 >= 1.6) Q2xBbin = 5;
			else if(x_B >= 0.28 && Q2 < 2.1) Q2xBbin = 6;
			else if(x_B >= 0.28 && Q2 >= 2.1) Q2xBbin = 7;

			IndexedList<double[]> tbin_l = new IndexedList<>(1);
			
			tbin_l.add(new double[] {0.23, 0.34, 0.51, 0.82}, 0);
			tbin_l.add(new double[] {0.27, 0.42, 0.59, 0.85}, 1);
			tbin_l.add(new double[] {0.25, 0.39, 0.56, 0.88}, 2);
			tbin_l.add(new double[] {0.35, 0.50, 0.63, 0.82}, 3);
			tbin_l.add(new double[] {0.34, 0.52, 0.70, 0.95}, 4);
			tbin_l.add(new double[] {0.31, 0.48, 0.68, 1.03}, 5);
			tbin_l.add(new double[] {0.46, 0.66, 0.85, 1.12}, 6);
			tbin_l.add(new double[] {0.56, 0.84, 1.14, 1.60}, 7);
			
			int ntbin = -1;
			if(-t_cal.mass2() < tbin_l.getItem(Q2xBbin)[0]) ntbin = 0;
			else if(-t_cal.mass2() >= tbin_l.getItem(Q2xBbin)[0] && -t_cal.mass2() < tbin_l.getItem(Q2xBbin)[1]) ntbin = 1;
			else if(-t_cal.mass2() >= tbin_l.getItem(Q2xBbin)[1] && -t_cal.mass2() < tbin_l.getItem(Q2xBbin)[2]) ntbin = 2;
			else if(-t_cal.mass2() >= tbin_l.getItem(Q2xBbin)[2] && -t_cal.mass2() < tbin_l.getItem(Q2xBbin)[3]) ntbin = 3;
			else if(-t_cal.mass2() >= tbin_l.getItem(Q2xBbin)[3]) ntbin = 4;
			
			if(ecount == 1 && protoncount == 1 && gammacount == 1
					&& v_e.pz() > -12 && v_e.pz() < 7 && Math.abs(v_e.pz()-v_proton.pz()) < 2.5 + (2.5/p_proton.p())
					&& fid_ecal_e == true
					&& (fid_ecal_gamma1 == true || fid_ftcal_gamma1 == true)
					&& p_gamma1.vect().theta(p_e.vect()) > 5
					&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30))
					&& 57.3*p_proton.theta() < 75)
			{
				if((det_gamma1 == 7
						&& p_gamma1.vect().theta(X_gamma1.vect()) < 4.5
						&& Math.sqrt((X.px()*X.px())+(X.py()*X.py())) < 0.3
						&& X.e() < 1.25
						&& X_proton.mass() > 0.4 && X_proton.mass() < 1.8)
						&& Q2 > 1 && W > 2
						&& -t_cal.mass2() > 0.938272*0.938272*x_B*x_B/(1-x_B+x_B*0.938272*0.938272/Q2*Q2))
				{
					{
						if(helic == 1)
						{
							hphi_neghel.fill(phi);
							histGroups_phi_neghel_tbin.getItem(Q2xBbin, ntbin).fill(phi);
						}
						if(helic == -1)
						{	
							hphi_poshel.fill(phi);
							histGroups_phi_poshel_tbin.getItem(Q2xBbin, ntbin).fill(phi);
						}
					}
				}
			}			
		}
	}
	
	static void processEvent_pi0(DataEvent event, int eventCounter)
	{
		if(event.hasBank("REC::Particle") && event.hasBank("REC::Track"))
		{
			HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
			int ecount = 0;
			int protoncount = 0;
			int gammacount = 0;
			int e_pindex = -1;
			boolean fid_ecal_e = false;
			int proton_pindex = -1;
			int gamma1_pindex = -1;
			int gamma2_pindex = -1;
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0.000511);
			LorentzVector p_proton = new LorentzVector(0, 0, 0, 0.938272);
			LorentzVector p_gamma1 = new LorentzVector(0, 0, 0, 0);
			LorentzVector p_gamma2 = new LorentzVector(0, 0, 0, 0);
			double E = 6.535;
			byte det_gamma1 = 0;
			byte det_gamma2 = 0;
			LorentzVector v_e = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_proton = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_gamma1 = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_gamma2 = new LorentzVector(0, 0, 0, 0);
			float cdchi2_proton = 100;
			short cdNDF_proton = 100;
			double cX = 0;
			double cY = 0;
			boolean fid_ecal_gamma1 = false;
			boolean fid_ecal_gamma2 = false;
			boolean fid_ftcal_gamma1 = false;
			boolean fid_ftcal_gamma2 = false;
			double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
			double[] cX_split  = new double[] {87, 82, 85, 77, 78, 82};
			double[] t_left  = new double[] {58.7356, 62.8204, 62.2296, 53.7756, 58.2888, 54.5822};
			double[] t_right  = new double[] {58.7477, 51.2589, 59.2357, 56.2415, 60.8219, 49.8914};
			double[] s_left  = new double[] {0.582053, 0.544976, 0.549788, 0.56899, 0.56414, 0.57343};
			double[] s_right  = new double[] {-0.591876, -0.562926, -0.562246, -0.563726, -0.568902, -0.550729};
			double[] r_left  = new double[] {64.9348, 64.7541, 67.832, 55.9324, 55.9225, 60.0997};
			double[] r_right  = new double[] {65.424, 54.6992, 63.6628, 57.8931, 56.5367, 56.4641};
			double[] q_left  = new double[] {0.745578, 0.606081, 0.729202, 0.627239, 0.503674, 0.717899};
			double[] q_right  = new double[] {-0.775022, -0.633863, -0.678901, -0.612458, -0.455319, -0.692481};
			for(int i = 0; i < rec.rows(); i++)
			{
				int pid = rec.getInt("pid", i);
				float px = rec.getFloat("px", i);
				float py = rec.getFloat("py", i);
				float pz = rec.getFloat("pz", i);
				float vx = rec.getFloat("vx", i);
				float vy = rec.getFloat("vy", i);
				float vz = rec.getFloat("vz", i);
				byte charge = rec.getByte("charge", i);
				float beta = rec.getFloat("beta", i);
				float p = (float) Math.sqrt((px*px)+(py*py)+(pz*pz));
				if(pid == 11)
				{
					ecount++;
					if(ecount == 1)
					{
						p_e = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.000511*0.000511)));
						v_e = new LorentzVector(vx, vy, vz, 0);
						e_pindex = i;
					}
				}
				if(pid == 2212)
				{
					protoncount++;
					if(protoncount == 1)
					{
						p_proton = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.938272*0.938272)));
						v_proton = new LorentzVector(vx, vy, vz, 0);
						proton_pindex = i;
					}
				}
				if(pid == 22)
				{
					gammacount++;
					if(gammacount == 1)
					{	
						p_gamma1 = new LorentzVector(px, py, pz, p);
						v_gamma1 = new LorentzVector(vx, vy, vz, 0);
						gamma1_pindex = i;
					}
					if(gammacount == 2)
					{
						p_gamma2 = new LorentzVector(px, py, pz, p);
						v_gamma2 = new LorentzVector(vx, vy, vz, 0);
						gamma2_pindex = i;
					}
				}
			}
			//Start: cal fiducial cut
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
					if(pindex == e_pindex && detector == 7 && layer == 1)
					{
						cX = x*Math.cos(rc*phi_rot[sector-1])+y*Math.sin(rc*phi_rot[sector-1]);
						cY = y*Math.cos(rc*phi_rot[sector-1])-x*Math.sin(rc*phi_rot[sector-1]);
						if(cX > cX_split[sector-1] && cY < s_left[sector-1]*(cX-t_left[sector-1])
							&& cY > s_right[sector-1]*(cX-t_right[sector-1]))
							fid_ecal_e = true;
						else if(cX < cX_split[sector-1] && cY < q_left[sector-1]*(cX-r_left[sector-1])
									&& cY > q_right[sector-1]*(cX-r_right[sector-1]))
							fid_ecal_e = true;
					}
					
					if(pindex == gamma1_pindex && detector == 7 && layer == 1)
					{
						cX = x*Math.cos(rc*phi_rot[sector-1])+y*Math.sin(rc*phi_rot[sector-1]);
						cY = y*Math.cos(rc*phi_rot[sector-1])-x*Math.sin(rc*phi_rot[sector-1]);
						if(cX > cX_split[sector-1] && cY < s_left[sector-1]*(cX-t_left[sector-1])
							&& cY > s_right[sector-1]*(cX-t_right[sector-1]))
							fid_ecal_gamma1 = true;
						else if(cX < cX_split[sector-1] && cY < q_left[sector-1]*(cX-r_left[sector-1])
									&& cY > q_right[sector-1]*(cX-r_right[sector-1]))
							fid_ecal_gamma1 = true;
						det_gamma1 = 7;
					}
					
					if(pindex == gamma1_pindex && detector == 7 && sector == 1)
					{
						if(layer == 1)
						{
							if((lw > 74.5 && lw < 79.5) || (lw > 85.5 && lw < 90.5) || (lw > 213 && lw < 218)
									|| (lw > 224.5 && lw < 229.5) || (lw > 410)) fid_ecal_gamma1 = false;
						}
						else if(layer == 4)
						{
							if((lv > 70 && lv < 93)) fid_ecal_gamma1 = false;
						}
						else if(layer == 7)
						{
							if((lu > 410.5)) fid_ecal_gamma1 = false;
							if((lv < 24) || (lv > 36.5 && lv < 60) || (lv > 411)) fid_ecal_gamma1 = false;
							if((lw < 21.5)) fid_ecal_gamma1 = false;
						}
					}
					else if(pindex == gamma1_pindex && detector == 7 && sector == 2)
					{
						if(layer == 1)
						{
							if((lv > 102 && lv < 113)) fid_ecal_gamma1 = false;
						}
						else if(layer == 4)
						{
							if((lu > 396)) fid_ecal_gamma1 = false;
							if((lw > 363)) fid_ecal_gamma1 = false;
						}
						else if(layer == 7)
						{
							if((lu < 12)) fid_ecal_gamma1 = false;
							if((lw < 10.5) || (lw > 376)) fid_ecal_gamma1 = false;
						}
					}
					else if(pindex == gamma1_pindex && detector == 7 && sector == 3)
					{
						if(layer == 1)
						{
							if((lv > 354.5 && lv < 376.5)) fid_ecal_gamma1 = false;
						}
						else if(layer == 4)
						{
							if((lu < 23)) fid_ecal_gamma1 = false;
							if((lw < 10) || (lw > 363)) fid_ecal_gamma1 = false;
						}
						else if(layer == 7)
						{
							if((lw > 387)) fid_ecal_gamma1 = false;
						}
					}
					else if(pindex == gamma1_pindex && detector == 7 && sector == 4)
					{
						if(layer == 1)
						{
							if((lv < 13) || (lv > 230 && lv < 240.5)) fid_ecal_gamma1 = false;
							if((lw > 410)) fid_ecal_gamma1 = false;
						}
						else if(layer == 4)
						{
							if((lu < 20.5)) fid_ecal_gamma1 = false;
						}
						else if(layer == 7)
						{
							if((lw < 32.5)) fid_ecal_gamma1 = false;
						}
					}
					else if(pindex == gamma1_pindex && detector == 7 && sector == 5)
					{
						if(layer == 4)
						{
							if((lv < 23)) fid_ecal_gamma1 = false;
							if((lw < 10)) fid_ecal_gamma1 = false;
						}
						else if(layer == 7)
						{
							if((lu > 193.5 && lu < 217)) fid_ecal_gamma1 = false;
							if((lv < 24)) fid_ecal_gamma1 = false;
						}
					}
					else if(pindex == gamma1_pindex && detector == 7 && sector == 6)
					{
						if(layer == 1)
						{
							if((lw > 174.5 && lw < 179.5) || (lw > 185.5 && lw < 190.5)) fid_ecal_gamma1 = false;
						}
						else if(layer == 4)
						{
							if((lv < 11.5)) fid_ecal_gamma1 = false;
							if((lu < 20.5)) fid_ecal_gamma1 = false;
						}
						else if(layer == 7)
						{
							if((lv < 12) || (lv > 423)) fid_ecal_gamma1 = false;
							if((lw < 32.5)) fid_ecal_gamma1 = false;
						}
					}

					if(pindex == gamma2_pindex && detector == 7 && layer == 1)
					{
						cX = x*Math.cos(rc*phi_rot[sector-1])+y*Math.sin(rc*phi_rot[sector-1]);
						cY = y*Math.cos(rc*phi_rot[sector-1])-x*Math.sin(rc*phi_rot[sector-1]);
						if(cX > cX_split[sector-1] && cY < s_left[sector-1]*(cX-t_left[sector-1])
							&& cY > s_right[sector-1]*(cX-t_right[sector-1]))
							fid_ecal_gamma2 = true;
						else if(cX < cX_split[sector-1] && cY < q_left[sector-1]*(cX-r_left[sector-1])
									&& cY > q_right[sector-1]*(cX-r_right[sector-1]))
							fid_ecal_gamma2 = true;
							det_gamma2 = 7;
					}
					
					if(pindex == gamma2_pindex && detector == 7 && sector == 1)
					{
						if(layer == 1)
						{
							if((lw > 74.5 && lw < 79.5) || (lw > 85.5 && lw < 90.5) || (lw > 213 && lw < 218)
									|| (lw > 224.5 && lw < 229.5) || (lw > 410)) fid_ecal_gamma2 = false;
						}
						else if(layer == 4)
						{
							if((lv > 70 && lv < 93)) fid_ecal_gamma2 = false;
						}
						else if(layer == 7)
						{
							if((lu > 410.5)) fid_ecal_gamma2 = false;
							if((lv < 24) || (lv > 36.5 && lv < 60) || (lv > 411)) fid_ecal_gamma2 = false;
							if((lw < 21.5)) fid_ecal_gamma2 = false;
						}
					}
					else if(pindex == gamma2_pindex && detector == 7 && sector == 2)
					{
						if(layer == 1)
						{
							if((lv > 102 && lv < 113)) fid_ecal_gamma2 = false;
						}
						else if(layer == 4)
						{
							if((lu > 396)) fid_ecal_gamma2 = false;
							if((lw > 363)) fid_ecal_gamma2 = false;
						}
						else if(layer == 7)
						{
							if((lu < 12)) fid_ecal_gamma2 = false;
							if((lw < 10.5) || (lw > 376)) fid_ecal_gamma2 = false;
						}
					}
					else if(pindex == gamma2_pindex && detector == 7 && sector == 3)
					{
						if(layer == 1)
						{
							if((lv > 354.5 && lv < 376.5)) fid_ecal_gamma2 = false;
						}
						else if(layer == 4)
						{
							if((lu < 23)) fid_ecal_gamma2 = false;
							if((lw < 10) || (lw > 363)) fid_ecal_gamma2 = false;
						}
						else if(layer == 7)
						{
							if((lw > 387)) fid_ecal_gamma2 = false;
						}
					}
					else if(pindex == gamma2_pindex && detector == 7 && sector == 4)
					{
						if(layer == 1)
						{
							if((lv < 13) || (lv > 230 && lv < 240.5)) fid_ecal_gamma2 = false;
							if((lw > 410)) fid_ecal_gamma2 = false;
						}
						else if(layer == 4)
						{
							if((lu < 20.5)) fid_ecal_gamma2 = false;
						}
						else if(layer == 7)
						{
							if((lw < 32.5)) fid_ecal_gamma2 = false;
						}
					}
					else if(pindex == gamma2_pindex && detector == 7 && sector == 5)
					{
						if(layer == 4)
						{
							if((lv < 23)) fid_ecal_gamma2 = false;
							if((lw < 10)) fid_ecal_gamma2 = false;
						}
						else if(layer == 7)
						{
							if((lu > 193.5 && lu < 217)) fid_ecal_gamma2 = false;
							if((lv < 24)) fid_ecal_gamma2 = false;
						}
					}
					else if(pindex == gamma2_pindex && detector == 7 && sector == 6)
					{
						if(layer == 1)
						{
							if((lw > 174.5 && lw < 179.5) || (lw > 185.5 && lw < 190.5)) fid_ecal_gamma2 = false;
						}
						else if(layer == 4)
						{
							if((lv < 11.5)) fid_ecal_gamma2 = false;
							if((lu < 20.5)) fid_ecal_gamma2 = false;
						}
						else if(layer == 7)
						{
							if((lv < 12) || (lv > 423)) fid_ecal_gamma2 = false;
							if((lw < 32.5)) fid_ecal_gamma2 = false;
						}
					}		
				}
			}
			
			if(event.hasBank("REC::ForwardTagger"))
			{
				HipoDataBank recft = (HipoDataBank) event.getBank("REC::ForwardTagger");
				for(int m = 0; m < recft.rows(); m++)
				{
					short pindex = recft.getShort("pindex", m);
					byte detector = recft.getByte("detector", m);
					byte layer = recft.getByte("layer", m);
					float x = recft.getFloat("x", m);
					float y = recft.getFloat("y", m);
					float z = recft.getFloat("z", m);
					Vector3D ftcal_det = new Vector3D (x, y, z);
					double xtrans = x + 18.36;
					double ytrans = y + 18.36;
					int cx = (int) (xtrans/1.53);
					int cy = (int) (ytrans/1.53);
					if(pindex == gamma1_pindex && detector == 10 && layer == 1)
					{
						if(ftcal_det.rho() > 8.8 && !(ftcal_det.rho() > 15.5)
							&& (x+8.5)*(x+8.5)+(y-10)*(y-10) > 1.5*1.5 
							&& (x+6)*(x+6)/(2.3*2.3)+(y+13)*(y+13)/(1.6*1.6) > 1
							&& (x-3.8)*(x-3.8)+(y+6.8)*(y+6.8) > 1.7*1.7) fid_ftcal_gamma1 = true;
						if(cx == 7 && cy == 15) fid_ftcal_gamma1 = false;
						det_gamma1 = 10;
					}
					if(pindex == gamma2_pindex && detector == 10 && layer == 1)
					{
						if(ftcal_det.rho() > 8.8 && !(ftcal_det.rho() > 15.5)
							&& (x+8.5)*(x+8.5)+(y-10)*(y-10) > 1.5*1.5 
							&& (x+6)*(x+6)/(2.3*2.3)+(y+13)*(y+13)/(1.6*1.6) > 1
							&& (x-3.8)*(x-3.8)+(y+6.8)*(y+6.8) > 1.7*1.7) fid_ftcal_gamma2 = true;
						if(cx == 7 && cy == 15) fid_ftcal_gamma2 = false;
						det_gamma2 = 10;
					}
				}
			}
			// End: cal fiducial cut
			HipoDataBank rectrac = (HipoDataBank) event.getBank("REC::Track");
			byte e_detector = 0;
			byte e_sector = 0;
			byte proton_detector = 0;
			byte gamma_sector = 0;
			byte X0_sector = 0;
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
			HipoDataBank recev = (HipoDataBank) event.getBank("REC::Event");
			byte helic = recev.getByte("helicity", 0);		
			LorentzVector q = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			LorentzVector qr = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			q.sub(p_e);
			qr.sub(p_e);
			double Q2 = -1*q.mass2();
			LorentzVector X = new LorentzVector(0, 0, 0, 0.938272);
			double x_B = Q2/(2*X.e()*q.e());
			LorentzVector t = new LorentzVector(p_proton.px(), p_proton.py(), p_proton.pz(), p_proton.e());
			t.sub(X);
			X.add(q);
			q.sub(p_gamma1);
			double W = X.mass();
			X.sub(p_proton);
			X.sub(p_gamma1);
			X.sub(p_gamma2);
			LorentzVector X_e = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_e.sub(p_proton);
			X_e.sub(p_gamma1);
			X_e.sub(p_gamma2);
			LorentzVector X_proton = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_proton.sub(p_e);
			X_proton.sub(p_gamma1);
			X_proton.sub(p_gamma2);
			LorentzVector X_tar = new LorentzVector(0, 0, 0, 0.938272);
			LorentzVector t_cal = new LorentzVector(X_proton.px(), X_proton.py(), X_proton.pz(), X_proton.e());
			t_cal.sub(X_tar);
			LorentzVector X_pi0 = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_pi0.sub(p_proton);
			X_pi0.sub(p_e);
			LorentzVector p_pi0 = new LorentzVector (0, 0, 0, 0);
			p_pi0.copy(p_gamma1);
			p_pi0.add(p_gamma2);
			LorentzVector beam = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			
			double trento = ((beam.vect().cross(p_e.vect())).dot(p_proton.vect()))/Math.abs((beam.vect().cross(p_e.vect())).dot(p_proton.vect()));
			double phi = trento*57.3*Math.acos(((beam.vect().cross(p_e.vect())).dot((p_proton.vect().cross(qr.vect()))))
							/((beam.vect().cross(p_e.vect()).mag()*(p_proton.vect().cross(qr.vect())).mag())));
			if(phi < 0) phi = phi+360;
			
			int Q2xBbin = -1;
			if(x_B < 0.16) Q2xBbin = 0;
			else if(x_B >= 0.16 && x_B < 0.21 && Q2 < 1.2) Q2xBbin = 1;
			else if(x_B >= 0.16 && x_B < 0.21 && Q2 >= 1.2) Q2xBbin = 2;
			else if(x_B >= 0.21 && x_B < 0.28 && Q2 < 1.2) Q2xBbin = 3;
			else if(x_B >= 0.21 && x_B < 0.28 && Q2 >= 1.2 && Q2 < 1.6) Q2xBbin = 4;
			else if(x_B >= 0.21 && x_B < 0.28 && Q2 >= 1.6) Q2xBbin = 5;
			else if(x_B >= 0.28 && Q2 < 2.1) Q2xBbin = 6;
			else if(x_B >= 0.28 && Q2 >= 2.1) Q2xBbin = 7;

			IndexedList<double[]> tbin_l = new IndexedList<>(1);
			
			tbin_l.add(new double[] {0.23, 0.34, 0.51, 0.82}, 0);
			tbin_l.add(new double[] {0.27, 0.42, 0.59, 0.85}, 1);
			tbin_l.add(new double[] {0.25, 0.39, 0.56, 0.88}, 2);
			tbin_l.add(new double[] {0.35, 0.50, 0.63, 0.82}, 3);
			tbin_l.add(new double[] {0.34, 0.52, 0.70, 0.95}, 4);
			tbin_l.add(new double[] {0.31, 0.48, 0.68, 1.03}, 5);
			tbin_l.add(new double[] {0.46, 0.66, 0.85, 1.12}, 6);
			tbin_l.add(new double[] {0.56, 0.84, 1.14, 1.60}, 7);
			
			int ntbin = -1;
			if(-t_cal.mass2() < tbin_l.getItem(Q2xBbin)[0]) ntbin = 0;
			else if(-t_cal.mass2() >= tbin_l.getItem(Q2xBbin)[0] && -t_cal.mass2() < tbin_l.getItem(Q2xBbin)[1]) ntbin = 1;
			else if(-t_cal.mass2() >= tbin_l.getItem(Q2xBbin)[1] && -t_cal.mass2() < tbin_l.getItem(Q2xBbin)[2]) ntbin = 2;
			else if(-t_cal.mass2() >= tbin_l.getItem(Q2xBbin)[2] && -t_cal.mass2() < tbin_l.getItem(Q2xBbin)[3]) ntbin = 3;
			else if(-t_cal.mass2() >= tbin_l.getItem(Q2xBbin)[3]) ntbin = 4;
			
			if(ecount == 1 && protoncount == 1 && gammacount == 2
					&& v_e.pz() > -12 && v_e.pz() < 7 && Math.abs(v_e.pz()-v_proton.pz()) < 2.5 + (2.5/p_proton.p())
					&& fid_ecal_e == true
					&& (fid_ecal_gamma1 == true || fid_ftcal_gamma1 == true)
					&& (fid_ecal_gamma2 == true || fid_ftcal_gamma2 == true)
					&& p_gamma1.vect().theta(p_e.vect()) > 5
					&& p_gamma2.vect().theta(p_e.vect()) > 5
					&& p_gamma1.vect().theta(p_gamma2.vect()) > 2
					&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30))
					&& 57.3*p_proton.theta() < 75)
			{
				if((det_gamma1 == 7 && det_gamma2 == 7
						&& p_pi0.vect().theta(X_pi0.vect())< 5.7
						&& Math.sqrt((X.px()*X.px())+(X.py()*X.py())) < 0.38
						&& X.e() < 1.25
						&& X_proton.mass() > 0.3 && X_proton.mass() < 1.9)
						&& Q2 > 1 && W > 2)
				{
					if(p_pi0.mass() > 0.103 && p_pi0.mass() < 0.163)
					{
						if(p_gamma1.vect().theta(p_gamma2.vect()) > 2.5)
						{
							if(helic == 1)
							{
								hphi_neghel_pi0.fill(phi);
								histGroups_phi_neghel_tbin_pi0.getItem(Q2xBbin, ntbin ).fill(phi);
							}
							if(helic == -1)
							{	
								hphi_poshel_pi0.fill(phi);
								histGroups_phi_poshel_tbin_pi0.getItem(Q2xBbin, ntbin ).fill(phi);
							}
						}
					}
				}
			}
		}
	}
	
	static void processEvent_pi0_sim(DataEvent event, int eventCounter)
	{
		if(event.hasBank("REC::Particle") && event.hasBank("REC::Track"))
		{
			HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
			int ecount = 0;
			int protoncount = 0;
			int gammacount = 0;
			int e_pindex = -1;
			boolean fid_ecal_e = false;
			int proton_pindex = -1;
			int gamma1_pindex = -1;
			int gamma2_pindex = -1;
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0.000511);
			LorentzVector p_proton = new LorentzVector(0, 0, 0, 0.938272);
			LorentzVector p_gamma1 = new LorentzVector(0, 0, 0, 0);
			LorentzVector p_gamma2 = new LorentzVector(0, 0, 0, 0);
			double E = 6.535;
			byte det_gamma1 = 0;
			float ftcalrad_gamma1 = 0;
			byte det_gamma2 = 0;
			float ftcalrad_gamma2 = 0;
			LorentzVector v_e = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_proton = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_gamma1 = new LorentzVector(0, 0, 0, 0);
			LorentzVector v_gamma2 = new LorentzVector(0, 0, 0, 0);
			float cdchi2_proton = 100;
			short cdNDF_proton = 100;
			double cX = 0;
			double cY = 0;
			boolean fid_ecal_gamma1 = false;
			boolean fid_ecal_gamma2 = false;
			boolean fid_ftcal_gamma1 = false;
			boolean fid_ftcal_gamma2 = false;
			double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
			double[] cX_split  = new double[] {87, 82, 85, 77, 78, 82};
			double[] t_left  = new double[] {58.7356, 62.8204, 62.2296, 53.7756, 58.2888, 54.5822};
			double[] t_right  = new double[] {58.7477, 51.2589, 59.2357, 56.2415, 60.8219, 49.8914};
			double[] s_left  = new double[] {0.582053, 0.544976, 0.549788, 0.56899, 0.56414, 0.57343};
			double[] s_right  = new double[] {-0.591876, -0.562926, -0.562246, -0.563726, -0.568902, -0.550729};
			double[] r_left  = new double[] {64.9348, 64.7541, 67.832, 55.9324, 55.9225, 60.0997};
			double[] r_right  = new double[] {65.424, 54.6992, 63.6628, 57.8931, 56.5367, 56.4641};
			double[] q_left  = new double[] {0.745578, 0.606081, 0.729202, 0.627239, 0.503674, 0.717899};
			double[] q_right  = new double[] {-0.775022, -0.633863, -0.678901, -0.612458, -0.455319, -0.692481};
			for(int i = 0; i < rec.rows(); i++)
			{
				int pid = rec.getInt("pid", i);
				float px = rec.getFloat("px", i);
				float py = rec.getFloat("py", i);
				float pz = rec.getFloat("pz", i);
				float vx = rec.getFloat("vx", i);
				float vy = rec.getFloat("vy", i);
				float vz = rec.getFloat("vz", i);
				byte charge = rec.getByte("charge", i);
				float beta = rec.getFloat("beta", i);
				float p = (float) Math.sqrt((px*px)+(py*py)+(pz*pz));
				if(pid == 11)
				{
					ecount++;
					if(ecount == 1)
					{
						p_e = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.000511*0.000511)));
						v_e = new LorentzVector(vx, vy, vz, 0);
						e_pindex = i;
					}
				}
				if(pid == 2212)
				{
					protoncount++;
					if(protoncount == 1)
					{
						p_proton = new LorentzVector(px, py, pz, Math.sqrt((p*p)+(0.938272*0.938272)));
						v_proton = new LorentzVector(vx, vy, vz, 0);
						proton_pindex = i;
					}
				}
				if(pid == 22)
				{
					gammacount++;
					if(gammacount == 1)
					{	
						p_gamma1 = new LorentzVector(px, py, pz, p);
						v_gamma1 = new LorentzVector(vx, vy, vz, 0);
						gamma1_pindex = i;
					}
					if(gammacount == 2)
					{
						p_gamma2 = new LorentzVector(px, py, pz, p);
						v_gamma2 = new LorentzVector(vx, vy, vz, 0);
						gamma2_pindex = i;
					}
				}
				//Start: cal fiducial cut
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
						if(pindex == e_pindex && detector == 7 && layer == 1)
						{
							cX = x*Math.cos(rc*phi_rot[sector-1])+y*Math.sin(rc*phi_rot[sector-1]);
							cY = y*Math.cos(rc*phi_rot[sector-1])-x*Math.sin(rc*phi_rot[sector-1]);
							if(cX > cX_split[sector-1] && cY < s_left[sector-1]*(cX-t_left[sector-1])
								&& cY > s_right[sector-1]*(cX-t_right[sector-1]))
								fid_ecal_e = true;
							else if(cX < cX_split[sector-1] && cY < q_left[sector-1]*(cX-r_left[sector-1])
										&& cY > q_right[sector-1]*(cX-r_right[sector-1]))
								fid_ecal_e = true;
						}
						
						if(pindex == gamma1_pindex && detector == 7 && layer == 1)
						{
							cX = x*Math.cos(rc*phi_rot[sector-1])+y*Math.sin(rc*phi_rot[sector-1]);
							cY = y*Math.cos(rc*phi_rot[sector-1])-x*Math.sin(rc*phi_rot[sector-1]);
							if(cX > cX_split[sector-1] && cY < s_left[sector-1]*(cX-t_left[sector-1])
								&& cY > s_right[sector-1]*(cX-t_right[sector-1]))
								fid_ecal_gamma1 = true;
							else if(cX < cX_split[sector-1] && cY < q_left[sector-1]*(cX-r_left[sector-1])
										&& cY > q_right[sector-1]*(cX-r_right[sector-1]))
								fid_ecal_gamma1 = true;
							det_gamma1 = 7;
						}
						
						if(pindex == gamma1_pindex && detector == 7 && sector == 1)
						{
							if(layer == 1)
							{
								if((lw > 74.5 && lw < 79.5) || (lw > 85.5 && lw < 90.5) || (lw > 213 && lw < 218)
										|| (lw > 224.5 && lw < 229.5) || (lw > 410)) fid_ecal_gamma1 = false;
							}
							else if(layer == 4)
							{
								if((lv > 70 && lv < 93)) fid_ecal_gamma1 = false;
							}
							else if(layer == 7)
							{
								if((lu > 410.5)) fid_ecal_gamma1 = false;
								if((lv < 24) || (lv > 36.5 && lv < 60) || (lv > 411)) fid_ecal_gamma1 = false;
								if((lw < 21.5)) fid_ecal_gamma1 = false;
							}
						}
						else if(pindex == gamma1_pindex && detector == 7 && sector == 2)
						{
							if(layer == 1)
							{
								if((lv > 102 && lv < 113)) fid_ecal_gamma1 = false;
							}
							else if(layer == 4)
							{
								if((lu > 396)) fid_ecal_gamma1 = false;
								if((lw > 363)) fid_ecal_gamma1 = false;
							}
							else if(layer == 7)
							{
								if((lu < 12)) fid_ecal_gamma1 = false;
								if((lw < 10.5) || (lw > 376)) fid_ecal_gamma1 = false;
							}
						}
						else if(pindex == gamma1_pindex && detector == 7 && sector == 3)
						{
							if(layer == 1)
							{
								if((lv > 354.5 && lv < 376.5)) fid_ecal_gamma1 = false;
							}
							else if(layer == 4)
							{
								if((lu < 23)) fid_ecal_gamma1 = false;
								if((lw < 10) || (lw > 363)) fid_ecal_gamma1 = false;
							}
							else if(layer == 7)
							{
								if((lw > 387)) fid_ecal_gamma1 = false;
							}
						}
						else if(pindex == gamma1_pindex && detector == 7 && sector == 4)
						{
							if(layer == 1)
							{
								if((lv < 13) || (lv > 230 && lv < 240.5)) fid_ecal_gamma1 = false;
								if((lw > 410)) fid_ecal_gamma1 = false;
							}
							else if(layer == 4)
							{
								if((lu < 20.5)) fid_ecal_gamma1 = false;
							}
							else if(layer == 7)
							{
								if((lw < 32.5)) fid_ecal_gamma1 = false;
							}
						}
						else if(pindex == gamma1_pindex && detector == 7 && sector == 5)
						{
							if(layer == 4)
							{
								if((lv < 23)) fid_ecal_gamma1 = false;
								if((lw < 10)) fid_ecal_gamma1 = false;
							}
							else if(layer == 7)
							{
								if((lu > 193.5 && lu < 217)) fid_ecal_gamma1 = false;
								if((lv < 24)) fid_ecal_gamma1 = false;
							}
						}
						else if(pindex == gamma1_pindex && detector == 7 && sector == 6)
						{
							if(layer == 1)
							{
								if((lw > 174.5 && lw < 179.5) || (lw > 185.5 && lw < 190.5)) fid_ecal_gamma1 = false;
							}
							else if(layer == 4)
							{
								if((lv < 11.5)) fid_ecal_gamma1 = false;
								if((lu < 20.5)) fid_ecal_gamma1 = false;
							}
							else if(layer == 7)
							{
								if((lv < 12) || (lv > 423)) fid_ecal_gamma1 = false;
								if((lw < 32.5)) fid_ecal_gamma1 = false;
							}
						}

						if(pindex == gamma2_pindex && detector == 7 && layer == 1)
						{
							cX = x*Math.cos(rc*phi_rot[sector-1])+y*Math.sin(rc*phi_rot[sector-1]);
							cY = y*Math.cos(rc*phi_rot[sector-1])-x*Math.sin(rc*phi_rot[sector-1]);
							if(cX > cX_split[sector-1] && cY < s_left[sector-1]*(cX-t_left[sector-1])
								&& cY > s_right[sector-1]*(cX-t_right[sector-1]))
								fid_ecal_gamma2 = true;
							else if(cX < cX_split[sector-1] && cY < q_left[sector-1]*(cX-r_left[sector-1])
										&& cY > q_right[sector-1]*(cX-r_right[sector-1]))
								fid_ecal_gamma2 = true;
								det_gamma2 = 7;
						}
						
						if(pindex == gamma2_pindex && detector == 7 && sector == 1)
						{
							if(layer == 1)
							{
								if((lw > 74.5 && lw < 79.5) || (lw > 85.5 && lw < 90.5) || (lw > 213 && lw < 218)
										|| (lw > 224.5 && lw < 229.5) || (lw > 410)) fid_ecal_gamma2 = false;
							}
							else if(layer == 4)
							{
								if((lv > 70 && lv < 93)) fid_ecal_gamma2 = false;
							}
							else if(layer == 7)
							{
								if((lu > 410.5)) fid_ecal_gamma2 = false;
								if((lv < 24) || (lv > 36.5 && lv < 60) || (lv > 411)) fid_ecal_gamma2 = false;
								if((lw < 21.5)) fid_ecal_gamma2 = false;
							}
						}
						else if(pindex == gamma2_pindex && detector == 7 && sector == 2)
						{
							if(layer == 1)
							{
								if((lv > 102 && lv < 113)) fid_ecal_gamma2 = false;
							}
							else if(layer == 4)
							{
								if((lu > 396)) fid_ecal_gamma2 = false;
								if((lw > 363)) fid_ecal_gamma2 = false;
							}
							else if(layer == 7)
							{
								if((lu < 12)) fid_ecal_gamma2 = false;
								if((lw < 10.5) || (lw > 376)) fid_ecal_gamma2 = false;
							}
						}
						else if(pindex == gamma2_pindex && detector == 7 && sector == 3)
						{
							if(layer == 1)
							{
								if((lv > 354.5 && lv < 376.5)) fid_ecal_gamma2 = false;
							}
							else if(layer == 4)
							{
								if((lu < 23)) fid_ecal_gamma2 = false;
								if((lw < 10) || (lw > 363)) fid_ecal_gamma2 = false;
							}
							else if(layer == 7)
							{
								if((lw > 387)) fid_ecal_gamma2 = false;
							}
						}
						else if(pindex == gamma2_pindex && detector == 7 && sector == 4)
						{
							if(layer == 1)
							{
								if((lv < 13) || (lv > 230 && lv < 240.5)) fid_ecal_gamma2 = false;
								if((lw > 410)) fid_ecal_gamma2 = false;
							}
							else if(layer == 4)
							{
								if((lu < 20.5)) fid_ecal_gamma2 = false;
							}
							else if(layer == 7)
							{
								if((lw < 32.5)) fid_ecal_gamma2 = false;
							}
						}
						else if(pindex == gamma2_pindex && detector == 7 && sector == 5)
						{
							if(layer == 4)
							{
								if((lv < 23)) fid_ecal_gamma2 = false;
								if((lw < 10)) fid_ecal_gamma2 = false;
							}
							else if(layer == 7)
							{
								if((lu > 193.5 && lu < 217)) fid_ecal_gamma2 = false;
								if((lv < 24)) fid_ecal_gamma2 = false;
							}
						}
						else if(pindex == gamma2_pindex && detector == 7 && sector == 6)
						{
							if(layer == 1)
							{
								if((lw > 174.5 && lw < 179.5) || (lw > 185.5 && lw < 190.5)) fid_ecal_gamma2 = false;
							}
							else if(layer == 4)
							{
								if((lv < 11.5)) fid_ecal_gamma2 = false;
								if((lu < 20.5)) fid_ecal_gamma2 = false;
							}
							else if(layer == 7)
							{
								if((lv < 12) || (lv > 423)) fid_ecal_gamma2 = false;
								if((lw < 32.5)) fid_ecal_gamma2 = false;
							}
						}		
					}
				}
				
				if(event.hasBank("REC::ForwardTagger"))
				{
					HipoDataBank recft = (HipoDataBank) event.getBank("REC::ForwardTagger");
					for(int m = 0; m < recft.rows(); m++)
					{
						short pindex = recft.getShort("pindex", m);
						byte detector = recft.getByte("detector", m);
						byte layer = recft.getByte("layer", m);
						float x = recft.getFloat("x", m);
						float y = recft.getFloat("y", m);
						float z = recft.getFloat("z", m);
						Vector3D ftcal_det = new Vector3D (x, y, z);
						double xtrans = x + 18.36;
						double ytrans = y + 18.36;
						int cx = (int) (xtrans/1.53);
						int cy = (int) (ytrans/1.53);
						if(pindex == gamma1_pindex && detector == 10 && layer == 1)
						{
							if(ftcal_det.rho() > 8.8 && !(ftcal_det.rho() > 15.5)
								&& (x+8.5)*(x+8.5)+(y-10)*(y-10) > 1.5*1.5 
								&& (x+6)*(x+6)/(2.3*2.3)+(y+13)*(y+13)/(1.6*1.6) > 1
								&& (x-3.8)*(x-3.8)+(y+6.8)*(y+6.8) > 1.7*1.7) fid_ftcal_gamma1 = true;
							if(cx == 7 && cy == 15) fid_ftcal_gamma1 = false;
							det_gamma1 = 10;
						}
						if(pindex == gamma2_pindex && detector == 10 && layer == 1)
						{
							if(ftcal_det.rho() > 8.8 && !(ftcal_det.rho() > 15.5)
								&& (x+8.5)*(x+8.5)+(y-10)*(y-10) > 1.5*1.5 
								&& (x+6)*(x+6)/(2.3*2.3)+(y+13)*(y+13)/(1.6*1.6) > 1
								&& (x-3.8)*(x-3.8)+(y+6.8)*(y+6.8) > 1.7*1.7) fid_ftcal_gamma2 = true;
							if(cx == 7 && cy == 15) fid_ftcal_gamma2 = false;
							det_gamma2 = 10;
						}
					}
				}
				// End: cal fiducial cut
			}
			HipoDataBank rectrac = (HipoDataBank) event.getBank("REC::Track");
			byte e_detector = 0;
			byte e_sector = 0;
			byte proton_detector = 0;
			byte gamma_sector = 0;
			byte X0_sector = 0;
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
			HipoDataBank recev = (HipoDataBank) event.getBank("REC::Event");
			byte helic = recev.getByte("helicity", 0);		
			LorentzVector q = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			LorentzVector qr = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			q.sub(p_e);
			qr.sub(p_e);
			double Q2 = -1*q.mass2();
			LorentzVector X_1gamma = new LorentzVector(0, 0, 0, 0.938272);
			double x_B = Q2/(2*X_1gamma.e()*q.e());
			LorentzVector t = new LorentzVector(p_proton.px(), p_proton.py(), p_proton.pz(), p_proton.e());
			t.sub(X_1gamma);
			X_1gamma.add(q);
			q.sub(p_gamma1);
			double W = X_1gamma.mass();
			X_1gamma.sub(p_proton);
			X_1gamma.sub(p_gamma1);
			LorentzVector X_2gamma = new LorentzVector(0, 0, 0, 0);
			X_2gamma.copy(X_1gamma);
			X_2gamma.sub(p_gamma2);
			LorentzVector X_e_1gamma = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_e_1gamma.sub(p_proton);
			X_e_1gamma.sub(p_gamma1);
			LorentzVector X_e_2gamma = new LorentzVector(0, 0, 0, 0);
			X_e_2gamma.copy(X_e_1gamma);
			X_e_2gamma.sub(p_gamma2);
			LorentzVector X_proton_1gamma = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_proton_1gamma.sub(p_e);
			X_proton_1gamma.sub(p_gamma1);
			LorentzVector X_proton_2gamma = new LorentzVector(0, 0, 0, 0);
			X_proton_2gamma.copy(X_proton_1gamma);
			X_proton_2gamma.sub(p_gamma2);
			LorentzVector X_tar = new LorentzVector(0, 0, 0, 0.938272);
			LorentzVector t_cal_1gamma = new LorentzVector(X_proton_1gamma.px(), X_proton_1gamma.py(), X_proton_1gamma.pz(), X_proton_1gamma.e());
			t_cal_1gamma.sub(X_tar);
			LorentzVector t_cal_2gamma = new LorentzVector(X_proton_2gamma.px(), X_proton_2gamma.py(), X_proton_2gamma.pz(), X_proton_2gamma.e());
			t_cal_2gamma.sub(X_tar);
			LorentzVector X_gamma1 = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_gamma1.sub(p_proton);
			X_gamma1.sub(p_e);
			LorentzVector X_pi0 = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511))+0.938272);
			X_pi0.sub(p_proton);
			X_pi0.sub(p_e);
			LorentzVector p_pi0 = new LorentzVector (0, 0, 0, 0);
			p_pi0.copy(p_gamma1);
			p_pi0.add(p_gamma2);
			LorentzVector beam = new LorentzVector(0, 0, E, Math.sqrt((E*E)+(0.000511*0.000511)));
			
			double trento = ((beam.vect().cross(p_e.vect())).dot(p_proton.vect()))/Math.abs((beam.vect().cross(p_e.vect())).dot(p_proton.vect()));
			double phi = trento*57.3*Math.acos(((beam.vect().cross(p_e.vect())).dot((p_proton.vect().cross(qr.vect()))))
							/((beam.vect().cross(p_e.vect()).mag()*(p_proton.vect().cross(qr.vect())).mag())));
			if(phi < 0) phi = phi+360;
			
			int Q2xBbin = -1;
			if(x_B < 0.16) Q2xBbin = 0;
			else if(x_B >= 0.16 && x_B < 0.21 && Q2 < 1.2) Q2xBbin = 1;
			else if(x_B >= 0.16 && x_B < 0.21 && Q2 >= 1.2) Q2xBbin = 2;
			else if(x_B >= 0.21 && x_B < 0.28 && Q2 < 1.2) Q2xBbin = 3;
			else if(x_B >= 0.21 && x_B < 0.28 && Q2 >= 1.2 && Q2 < 1.6) Q2xBbin = 4;
			else if(x_B >= 0.21 && x_B < 0.28 && Q2 >= 1.6) Q2xBbin = 5;
			else if(x_B >= 0.28 && Q2 < 2.1) Q2xBbin = 6;
			else if(x_B >= 0.28 && Q2 >= 2.1) Q2xBbin = 7;

			IndexedList<double[]> tbin_l = new IndexedList<>(1);
			
			tbin_l.add(new double[] {0.23, 0.34, 0.51, 0.82}, 0);
			tbin_l.add(new double[] {0.27, 0.42, 0.59, 0.85}, 1);
			tbin_l.add(new double[] {0.25, 0.39, 0.56, 0.88}, 2);
			tbin_l.add(new double[] {0.35, 0.50, 0.63, 0.82}, 3);
			tbin_l.add(new double[] {0.34, 0.52, 0.70, 0.95}, 4);
			tbin_l.add(new double[] {0.31, 0.48, 0.68, 1.03}, 5);
			tbin_l.add(new double[] {0.46, 0.66, 0.85, 1.12}, 6);
			tbin_l.add(new double[] {0.56, 0.84, 1.14, 1.60}, 7);
			
			int ntbin_1gamma = -1;
			if(-t_cal_1gamma.mass2() < tbin_l.getItem(Q2xBbin)[0]) ntbin_1gamma = 0;
			else if(-t_cal_1gamma.mass2() >= tbin_l.getItem(Q2xBbin)[0] && -t_cal_1gamma.mass2() < tbin_l.getItem(Q2xBbin)[1]) ntbin_1gamma = 1;
			else if(-t_cal_1gamma.mass2() >= tbin_l.getItem(Q2xBbin)[1] && -t_cal_1gamma.mass2() < tbin_l.getItem(Q2xBbin)[2]) ntbin_1gamma = 2;
			else if(-t_cal_1gamma.mass2() >= tbin_l.getItem(Q2xBbin)[2] && -t_cal_1gamma.mass2() < tbin_l.getItem(Q2xBbin)[3]) ntbin_1gamma = 3;
			else if(-t_cal_1gamma.mass2() >= tbin_l.getItem(Q2xBbin)[3]) ntbin_1gamma = 4;
			
			int ntbin_2gamma = -1;
			if(-t_cal_2gamma.mass2() < tbin_l.getItem(Q2xBbin)[0]) ntbin_2gamma = 0;
			else if(-t_cal_2gamma.mass2() >= tbin_l.getItem(Q2xBbin)[0] && -t_cal_2gamma.mass2() < tbin_l.getItem(Q2xBbin)[1]) ntbin_2gamma = 1;
			else if(-t_cal_2gamma.mass2() >= tbin_l.getItem(Q2xBbin)[1] && -t_cal_2gamma.mass2() < tbin_l.getItem(Q2xBbin)[2]) ntbin_2gamma = 2;
			else if(-t_cal_2gamma.mass2() >= tbin_l.getItem(Q2xBbin)[2] && -t_cal_2gamma.mass2() < tbin_l.getItem(Q2xBbin)[3]) ntbin_2gamma = 3;
			else if(-t_cal_2gamma.mass2() >= tbin_l.getItem(Q2xBbin)[3]) ntbin_2gamma = 4;
			
			if(ecount == 1 && protoncount == 1 && gammacount == 1
					&& v_e.pz() > -12 && v_e.pz() < 7 && Math.abs(v_e.pz()-v_proton.pz()) < 2.5 + (2.5/p_proton.p())
					&& fid_ecal_e == true
					&& (fid_ecal_gamma1 == true || fid_ftcal_gamma1 == true)
					&& p_gamma1.vect().theta(p_e.vect()) > 5
					&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30))
					&& 57.3*p_proton.theta() < 75)
			{
				if((det_gamma1 == 7
						&& p_gamma1.vect().theta(X_gamma1.vect()) < 4.5
						&& Math.sqrt((X_1gamma.px()*X_1gamma.px())+(X_1gamma.py()*X_1gamma.py())) < 0.3
						&& X_1gamma.e() < 1.25
						&& X_proton_1gamma.mass() > 0.4 && X_proton_1gamma.mass() < 1.8)
						&& Q2 > 1 && W > 2
						&& -t_cal_1gamma.mass2() > 0.938272*0.938272*x_B*x_B/(1-x_B+x_B*0.938272*0.938272/Q2*Q2))
				{
					{
						{
							hphi_1gamma_pi0_sim.fill(phi);
							histGroups_phi_1gamma_pi0_sim_tbin.getItem(Q2xBbin, ntbin_1gamma).fill(phi);
						}
					}
				}
			}
			if(ecount == 1 && protoncount == 1 && gammacount == 2
					&& v_e.pz() > -12 && v_e.pz() < 7 && Math.abs(v_e.pz()-v_proton.pz()) < 2.5 + (2.5/p_proton.p())
					&& fid_ecal_e == true
					&& (fid_ecal_gamma1 == true || fid_ftcal_gamma1 == true)
					&& (fid_ecal_gamma2 == true || fid_ftcal_gamma2 == true)
					&& p_gamma1.vect().theta(p_e.vect()) > 5
					&& p_gamma2.vect().theta(p_e.vect()) > 5
					&& p_gamma1.vect().theta(p_gamma2.vect()) > 2
					&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30))
					&& 57.3*p_proton.theta() < 75)
			{
				if((det_gamma1 == 7 && det_gamma2 == 7
						&& p_pi0.vect().theta(X_pi0.vect())< 5.7
						&& Math.sqrt((X_2gamma.px()*X_2gamma.px())+(X_2gamma.py()*X_2gamma.py())) < 0.38
						&& X_2gamma.e() < 1.25
						&& X_proton_2gamma.mass() > 0.3 && X_proton_2gamma.mass() < 1.9)
						&& Q2 > 1 && W > 2)
				{
					if(p_pi0.mass() > 0.103 && p_pi0.mass() < 0.163)
					{
						if(p_gamma1.vect().theta(p_gamma2.vect()) > 2.5)
						{
							{
								hphi_2gamma_pi0_sim.fill(phi);
								histGroups_phi_2gamma_pi0_sim_tbin.getItem(Q2xBbin, ntbin_2gamma).fill(phi);
							}
						}
					}
				}
			}
		}
	}

	public static void main(String[] args) {
		
		phi_poshel_tbin_histos();
		phi_neghel_tbin_histos();
		phi_1gamma_pi0_sim_tbin_histos();
		phi_2gamma_pi0_sim_tbin_histos();
		
		reader.open("C:/Users/joshtanj/Documents/download/skim_epg_bank_merged_6535MeV_skim16.hipo");
		reader_pi0.open("C:/Users/joshtanj/Documents/download/skim_epgg_bank_merged_6535MeV_skim18.hipo");
		reader_pi0_sim.open("C:/Users/joshtanj/Documents/download/skim_epgg_epg_bank_rgk_fall2018_FTOff_noM_OUT.hipo");
		
		int eventCounter = 1;
		while(reader.hasEvent())// && eventCounter < 10000000)
		{
			processEvent(reader.getNextEvent(), eventCounter);
			if(eventCounter%50000 == 0) System.out.println("Event: " + eventCounter);
			eventCounter++;
		}
		int Nevent = eventCounter-1;
		System.out.println("Number of events: " + Nevent);
		
		int eventCounter_pi0 = 1;
		while(reader_pi0.hasEvent())// && eventCounter_pi0 < 10000000)
		{
			processEvent_pi0(reader_pi0.getNextEvent(), eventCounter_pi0);
			if(eventCounter_pi0%50000 == 0) System.out.println("DV#pi^0P event: " + eventCounter_pi0);
			eventCounter_pi0++;
		}
		int Nevent_pi0 = eventCounter_pi0-1;
		System.out.println("Number of DV#pi^0P events: " + Nevent_pi0);
		
		int eventCounter_pi0_sim = 1;
		while(reader_pi0_sim.hasEvent())// && eventCounter_pi0_sim < 10000000)
		{
			processEvent_pi0_sim(reader_pi0_sim.getNextEvent(), eventCounter_pi0_sim);
			if(eventCounter_pi0_sim%50000 == 0) System.out.println("DV#pi^0P_sim event: " + eventCounter_pi0_sim);
			eventCounter_pi0_sim++;
		}
		int Nevent_pi0_sim = eventCounter_pi0_sim-1;
		System.out.println("Number of DV#pi^0P_sim events: " + Nevent_pi0_sim);
				
		List<Double> bsa = new ArrayList<>();
		List<Double> phi_rb_poshel = new ArrayList<>();
		List<Double> phi_rb_neghel = new ArrayList<>();
		List<Double> bsa_pi0 = new ArrayList<>();
		List<Double> phi_rb_poshel_pi0 = new ArrayList<>();
		List<Double> phi_rb_neghel_pi0 = new ArrayList<>();
		List<Double> phi_rb_1gamma_pi0_sim = new ArrayList<>();
		List<Double> phi_rb_2gamma_pi0_sim = new ArrayList<>();
		List<Double> r_pi0_sim = new ArrayList<>();
		List<Double> R_pi0 = new ArrayList<>();
		List<Double> f_pi0_sim = new ArrayList<>();
		List<Double> bsa_sub_pi0_sim = new ArrayList<>();
		
		for(int prb = 0; prb < 20; prb++)
		{
			phi_rb_poshel.add((double) 0);
			phi_rb_neghel.add((double) 0);
			phi_rb_poshel_pi0.add((double) 0);
			phi_rb_neghel_pi0.add((double) 0);
			phi_rb_1gamma_pi0_sim.add((double) 0);
			phi_rb_2gamma_pi0_sim.add((double) 0);
		}
		
		int[] pb = new int[] {0, 5, 10, 15, 20, 26, 32, 38, 44, 52, 60,
								68, 76, 82, 88, 94, 100, 105, 110, 115, 120};
		
		for(int rb = 0; rb < 20; rb++)
		{
			for(int rbent = pb[rb]; rbent < pb[rb+1]; rbent++)
			{
				double phi_poshel_temp = phi_rb_poshel.get(rb);
				phi_rb_poshel.set(rb, phi_poshel_temp + hphi_poshel.getBinContent(rbent));
				double phi_neghel_temp = phi_rb_neghel.get(rb);
				phi_rb_neghel.set(rb, phi_neghel_temp + hphi_neghel.getBinContent(rbent));
				double phi_poshel_temp_pi0 = phi_rb_poshel_pi0.get(rb);
				phi_rb_poshel_pi0.set(rb, phi_poshel_temp_pi0 + hphi_poshel_pi0.getBinContent(rbent));
				double phi_neghel_temp_pi0 = phi_rb_neghel_pi0.get(rb);
				phi_rb_neghel_pi0.set(rb, phi_neghel_temp_pi0 + hphi_neghel_pi0.getBinContent(rbent));
				double phi_1gamma_pi0_sim_temp = phi_rb_1gamma_pi0_sim.get(rb);
				phi_rb_1gamma_pi0_sim.set(rb, phi_1gamma_pi0_sim_temp + hphi_1gamma_pi0_sim.getBinContent(rbent));
				double phi_2gamma_pi0_sim_temp = phi_rb_2gamma_pi0_sim.get(rb);
				phi_rb_2gamma_pi0_sim.set(rb, phi_2gamma_pi0_sim_temp + hphi_2gamma_pi0_sim.getBinContent(rbent));				
			}
		}
		
		for(int i = 0; i < 20; i++)
		{
			bsa.add((1/0.863)*(phi_rb_poshel.get(i)-phi_rb_neghel.get(i))/(phi_rb_poshel.get(i)+phi_rb_neghel.get(i)));
			bsa_pi0.add((1/0.863)*(phi_rb_poshel_pi0.get(i)-phi_rb_neghel_pi0.get(i))/(phi_rb_poshel_pi0.get(i)+phi_rb_neghel_pi0.get(i)));
			r_pi0_sim.add(phi_rb_1gamma_pi0_sim.get(i)/phi_rb_2gamma_pi0_sim.get(i));
			R_pi0.add((phi_rb_poshel_pi0.get(i)+phi_rb_neghel_pi0.get(i))
					/(phi_rb_poshel.get(i)+phi_rb_neghel.get(i)));
			f_pi0_sim.add(R_pi0.get(i)*r_pi0_sim.get(i));
			bsa_sub_pi0_sim.add((bsa.get(i)-(f_pi0_sim.get(i)*bsa_pi0.get(i)))/(1-f_pi0_sim.get(i)));
		}
		
		List<Double> bsaerr = new ArrayList<>();
		List<Double> bsaerr_pi0 = new ArrayList<>();
		List<Double> rerr_pi0_sim = new ArrayList<>();
		List<Double> ferr_pi0_sim = new ArrayList<>();
		List<Double> bsaerr_sub_pi0_sim = new ArrayList<>();
		for(int j = 0; j < 20; j++)
		{
			bsaerr.add((1/0.863)*Math.sqrt((1-(0.863*0.863*bsa.get(j)*bsa.get(j)))/(phi_rb_poshel.get(j)+phi_rb_neghel.get(j))));
			bsaerr_pi0.add((1/0.863)*Math.sqrt((1-(0.863*0.863*bsa_pi0.get(j)*bsa_pi0.get(j)))/(phi_rb_poshel_pi0.get(j)+phi_rb_neghel_pi0.get(j))));
			rerr_pi0_sim.add(Math.sqrt(r_pi0_sim.get(j)*(1+r_pi0_sim.get(j))/phi_rb_2gamma_pi0_sim.get(j)));
			ferr_pi0_sim.add(r_pi0_sim.get(j)*Math.sqrt(R_pi0.get(j)*(1+R_pi0.get(j))
							/(phi_rb_poshel.get(j)+phi_rb_neghel.get(j)))+R_pi0.get(j)*rerr_pi0_sim.get(j));
			double lbsa = bsa.get(j);
			double lbsa_pi0 = bsa_pi0.get(j);
			double lbsaerr = bsaerr.get(j);
			double lbsaerr_pi0 = bsaerr_pi0.get(j);
			double lf = f_pi0_sim.get(j);
			double lferr = ferr_pi0_sim.get(j);
			bsaerr_sub_pi0_sim.add(Math.sqrt((1+lf*lf)*(lbsaerr*lbsaerr)+lf*lf*(1+lf*lf)*(lbsaerr_pi0*lbsaerr_pi0)
										+(lbsa*lbsa+lbsa_pi0*lbsa_pi0)*lferr*lferr)/((1-lf)*(1-lf)));
		}
		
		double[] phibc = new double[] {7.5, 22.5, 37.5, 52.5, 69, 87, 105, 123, 144, 168,
										192, 206, 227, 245, 263, 281, 297.5, 312.5, 337.5, 352.5};
		
		double pentmin = 0;
		double pentmin_pi0 = 0;

		GraphErrors hBSAy = new GraphErrors();
		GraphErrors hBSAy_pi0 = new GraphErrors();
		GraphErrors hf_pi0_sim = new GraphErrors();
		GraphErrors hBSAy_sub_pi0_sim = new GraphErrors();
		
		for(int k = 0; k < 20; k++)
		{	
			if(phi_rb_poshel.get(k)+phi_rb_neghel.get(k) > pentmin) hBSAy.addPoint(phibc[k], bsa.get(k), 0, bsaerr.get(k));
			if(phi_rb_poshel_pi0.get(k)+phi_rb_neghel_pi0.get(k) > pentmin_pi0) hBSAy_pi0.addPoint(phibc[k], bsa_pi0.get(k), 0, bsaerr_pi0.get(k));
			if(phi_rb_poshel.get(k)+phi_rb_neghel.get(k) > pentmin && phi_rb_poshel_pi0.get(k)+phi_rb_neghel_pi0.get(k) > pentmin_pi0)
			{
				hf_pi0_sim.addPoint(phibc[k], f_pi0_sim.get(k), 0, ferr_pi0_sim.get(k));
				hBSAy_sub_pi0_sim.addPoint(phibc[k], bsa_sub_pi0_sim.get(k), 0, bsaerr_sub_pi0_sim.get(k));
			}
		}
		
		IndexedList<Double> bsa_tbin = new IndexedList<>(3);
		IndexedList<Double> phi_rb_poshel_tbin = new IndexedList<>(3);
		IndexedList<Double> phi_rb_neghel_tbin = new IndexedList<>(3);
		IndexedList<Double> bsa_tbin_pi0 = new IndexedList<>(3);
		IndexedList<Double> phi_rb_poshel_tbin_pi0 = new IndexedList<>(3);
		IndexedList<Double> phi_rb_neghel_tbin_pi0 = new IndexedList<>(3);
		IndexedList<Double> phi_rb_1gamma_pi0_sim_tbin = new IndexedList<>(3);
		IndexedList<Double> phi_rb_2gamma_pi0_sim_tbin = new IndexedList<>(3);
		IndexedList<Double> r_pi0_sim_tbin = new IndexedList<>(3);
		IndexedList<Double> R_pi0_tbin = new IndexedList<>(3);
		IndexedList<Double> f_pi0_sim_tbin = new IndexedList<>(3);
		IndexedList<Double> bsa_sub_pi0_sim_tbin = new IndexedList<>(3);
		
		for(int prbbin = 0; prbbin  < 8; prbbin++)
		{
			for (int prbtbin = 0; prbtbin < 5; prbtbin++)
			{
				for(int prbent = 0; prbent < 20; prbent++)
				{
					phi_rb_poshel_tbin.add((double) 0, prbbin, prbtbin, prbent);
					phi_rb_neghel_tbin.add((double) 0, prbbin, prbtbin, prbent);
					phi_rb_poshel_tbin_pi0.add((double) 0, prbbin, prbtbin, prbent);
					phi_rb_neghel_tbin_pi0.add((double) 0, prbbin, prbtbin, prbent);
					phi_rb_1gamma_pi0_sim_tbin.add((double) 0, prbbin, prbtbin, prbent);
					phi_rb_2gamma_pi0_sim_tbin.add((double) 0, prbbin, prbtbin, prbent);
				}
			}
		}
		
		for(int pbin = 0; pbin  < 8; pbin++)
		{
			for (int ptbin = 0; ptbin < 5; ptbin++)
			{
				for(int rb = 0; rb < 20; rb++)
				{
					for(int rbent = pb[rb]; rbent < pb[rb+1]; rbent++)
					{
						double phi_poshel_temp = phi_rb_poshel_tbin.getItem(pbin, ptbin, rb);
						phi_rb_poshel_tbin.add(phi_poshel_temp + histGroups_phi_poshel_tbin.getItem(pbin, ptbin).getBinContent(rbent),
												pbin, ptbin, rb);
						double phi_neghel_temp = phi_rb_neghel_tbin.getItem(pbin, ptbin, rb);
						phi_rb_neghel_tbin.add(phi_neghel_temp + histGroups_phi_neghel_tbin.getItem(pbin, ptbin).getBinContent(rbent),
												pbin, ptbin, rb);
						double phi_poshel_temp_pi0 = phi_rb_poshel_tbin_pi0.getItem(pbin, ptbin, rb);
						phi_rb_poshel_tbin_pi0.add(phi_poshel_temp_pi0 + histGroups_phi_poshel_tbin_pi0.getItem(pbin, ptbin).getBinContent(rbent),
												pbin, ptbin, rb);
						double phi_neghel_temp_pi0 = phi_rb_neghel_tbin_pi0.getItem(pbin, ptbin, rb);
						phi_rb_neghel_tbin_pi0.add(phi_neghel_temp_pi0 + histGroups_phi_neghel_tbin_pi0.getItem(pbin, ptbin).getBinContent(rbent),
												pbin, ptbin, rb);
						double phi_1gamma_pi0_sim_temp = phi_rb_1gamma_pi0_sim_tbin.getItem(pbin, ptbin, rb);
						phi_rb_1gamma_pi0_sim_tbin.add(phi_1gamma_pi0_sim_temp + histGroups_phi_1gamma_pi0_sim_tbin.getItem(pbin, ptbin).getBinContent(rbent),
														pbin, ptbin, rb);
						double phi_2gamma_pi0_sim_temp = phi_rb_2gamma_pi0_sim_tbin.getItem(pbin, ptbin, rb);
						phi_rb_2gamma_pi0_sim_tbin.add(phi_2gamma_pi0_sim_temp + histGroups_phi_2gamma_pi0_sim_tbin.getItem(pbin, ptbin).getBinContent(rbent),
														pbin, ptbin, rb);
					}
				}
			}
		}
		
		for(int ibin = 0; ibin  < 8; ibin++)
		{
			for (int itbin = 0; itbin < 5; itbin++)
			{
				for(int ient = 0; ient < 20; ient++)
				{
					bsa_tbin.add((1/0.863)*(phi_rb_poshel_tbin.getItem(ibin, itbin, ient) - phi_rb_neghel_tbin.getItem(ibin, itbin, ient))
									/(phi_rb_poshel_tbin.getItem(ibin, itbin, ient) + phi_rb_neghel_tbin.getItem(ibin, itbin, ient)), ibin, itbin, ient);
					bsa_tbin_pi0.add((1/0.863)*(phi_rb_poshel_tbin_pi0.getItem(ibin, itbin, ient) - phi_rb_neghel_tbin_pi0.getItem(ibin, itbin, ient))
										/(phi_rb_poshel_tbin_pi0.getItem(ibin, itbin, ient) + phi_rb_neghel_tbin_pi0.getItem(ibin, itbin, ient)),
										ibin, itbin, ient);
					r_pi0_sim_tbin.add(phi_rb_1gamma_pi0_sim_tbin.getItem(ibin, itbin, ient)/phi_rb_2gamma_pi0_sim_tbin.getItem(ibin, itbin, ient),
										ibin, itbin, ient);
					R_pi0_tbin.add((phi_rb_poshel_tbin_pi0.getItem(ibin, itbin, ient) + phi_rb_neghel_tbin_pi0.getItem(ibin, itbin, ient))
										/(phi_rb_poshel_tbin.getItem(ibin, itbin, ient) + phi_rb_neghel_tbin.getItem(ibin, itbin, ient)),
										ibin, itbin, ient);
					f_pi0_sim_tbin.add(R_pi0_tbin.getItem(ibin, itbin, ient)*r_pi0_sim_tbin.getItem(ibin, itbin, ient), ibin, itbin, ient);
					bsa_sub_pi0_sim_tbin.add((bsa_tbin.getItem(ibin, itbin, ient)-(f_pi0_sim_tbin.getItem(ibin, itbin, ient)
												*bsa_tbin_pi0.getItem(ibin, itbin, ient)))/(1-f_pi0_sim_tbin.getItem(ibin, itbin, ient)), ibin, itbin, ient);
				}
			}
		}
		
		IndexedList<Double> bsaerr_tbin = new IndexedList<>(3);
		IndexedList<Double> bsaerr_tbin_pi0 = new IndexedList<>(3);
		IndexedList<Double> rerr_pi0_sim_tbin = new IndexedList<>(3);
		IndexedList<Double> ferr_pi0_sim_tbin = new IndexedList<>(3);
		IndexedList<Double> bsaerr_sub_pi0_sim_tbin = new IndexedList<>(3);
		
		for(int jbin = 0; jbin  < 8; jbin++)
		{
			for (int jtbin = 0; jtbin < 5; jtbin++)
			{
				for(int jent = 0; jent < 20; jent++)
				{
					bsaerr_tbin.add((1/0.863)*Math.sqrt((1-(0.863*0.863*bsa_tbin.getItem(jbin, jtbin, jent)*bsa_tbin.getItem(jbin, jtbin, jent)))
									/(phi_rb_poshel_tbin.getItem(jbin, jtbin, jent) + phi_rb_neghel_tbin.getItem(jbin, jtbin, jent))), jbin, jtbin, jent);
					bsaerr_tbin_pi0.add((1/0.863)*Math.sqrt((1-(0.863*0.863*bsa_tbin_pi0.getItem(jbin, jtbin, jent)*bsa_tbin_pi0.getItem(jbin, jtbin, jent)))
										/(phi_rb_poshel_tbin_pi0.getItem(jbin, jtbin, jent) + phi_rb_neghel_tbin_pi0.getItem(jbin, jtbin, jent))),
										jbin, jtbin, jent);
					rerr_pi0_sim_tbin.add(Math.sqrt(r_pi0_sim_tbin.getItem(jbin, jtbin, jent)*(1+r_pi0_sim_tbin.getItem(jbin, jtbin, jent))
											/phi_rb_2gamma_pi0_sim_tbin.getItem(jbin, jtbin, jent)), jbin, jtbin, jent);
					ferr_pi0_sim_tbin.add(r_pi0_sim_tbin.getItem(jbin, jtbin, jent)*Math.sqrt(R_pi0_tbin.getItem(jbin, jtbin, jent)
											*(1+R_pi0_tbin.getItem(jbin, jtbin, jent))/(phi_rb_poshel_tbin.getItem(jbin, jtbin, jent)
											+phi_rb_neghel_tbin.getItem(jbin, jtbin, jent)))+R_pi0_tbin.getItem(jbin, jtbin, jent)
											*rerr_pi0_sim_tbin.getItem(jbin, jtbin, jent), jbin, jtbin, jent);
					double lbsa = bsa_tbin.getItem(jbin, jtbin, jent);
					double lbsa_pi0 = bsa_tbin_pi0.getItem(jbin, jtbin, jent);
					double lbsaerr = bsaerr_tbin.getItem(jbin, jtbin, jent);
					double lbsaerr_pi0 = bsaerr_tbin_pi0.getItem(jbin, jtbin, jent);
					double lf = f_pi0_sim_tbin.getItem(jbin, jtbin, jent);
					double lferr = ferr_pi0_sim_tbin.getItem(jbin, jtbin, jent);
					bsaerr_sub_pi0_sim_tbin.add(Math.sqrt((1+lf*lf)*(lbsaerr*lbsaerr)+lf*lf*(1+lf*lf)*(lbsaerr_pi0*lbsaerr_pi0)
												+(lbsa*lbsa+lbsa_pi0*lbsa_pi0)*lferr*lferr)/((1-lf)*(1-lf)), jbin, jtbin, jent);
				}
			}
		}
				
		IndexedList<GraphErrors> hBSAy_tbin = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> halpha_vs_t = new IndexedList<GraphErrors>(1);
		IndexedList<GraphErrors> hBSAy_tbin_pi0 = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> halpha_vs_t_pi0 = new IndexedList<GraphErrors>(1);
		IndexedList<GraphErrors> hf_pi0_sim_tbin = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> hBSAy_sub_pi0_sim_tbin = new IndexedList<GraphErrors>(2);
		IndexedList<GraphErrors> halpha_corr_pi0_sim_vs_t = new IndexedList<GraphErrors>(1);
		
		IndexedList<double[]> tba = new IndexedList<>(1);
		
		tba.add(new double[] {0.16, 0.29, 0.42, 0.64, 1.25}, 0);
		tba.add(new double[] {0.19, 0.34, 0.50, 0.71, 1.19}, 1);
		tba.add(new double[] {0.17, 0.32, 0.47, 0.71, 1.33}, 2);
		tba.add(new double[] {0.24, 0.43, 0.57, 0.71, 1.08}, 3);
		tba.add(new double[] {0.23, 0.44, 0.61, 0.82, 1.30}, 4);
		tba.add(new double[] {0.21, 0.39, 0.58, 0.83, 1.49}, 5);
		tba.add(new double[] {0.33, 0.57, 0.75, 0.97, 1.50}, 6);
		tba.add(new double[] {0.39, 0.70, 0.98, 1.34, 2.07}, 7);
		
		for(int kbin = 0; kbin  < 8; kbin++)
		{
			GraphErrors halpha_vs_t_ini = new GraphErrors();
			halpha_vs_t.add(halpha_vs_t_ini, kbin);
			GraphErrors halpha_vs_t_ini_pi0 = new GraphErrors();
			halpha_vs_t_pi0.add(halpha_vs_t_ini_pi0, kbin);
			GraphErrors halpha_corr_pi0_sim_vs_t_ini = new GraphErrors();
			halpha_corr_pi0_sim_vs_t.add(halpha_corr_pi0_sim_vs_t_ini, kbin);
			for(int ktbin = 0; ktbin < 5; ktbin++)
			{
				GraphErrors hBSAy_tbin_ini = new GraphErrors();
				hBSAy_tbin.add(hBSAy_tbin_ini, kbin, ktbin);
				GraphErrors hBSAy_tbin_ini_pi0 = new GraphErrors();
				hBSAy_tbin_pi0.add(hBSAy_tbin_ini_pi0, kbin, ktbin);
				GraphErrors hf_pi0_sim_tbin_ini = new GraphErrors();
				hf_pi0_sim_tbin.add(hf_pi0_sim_tbin_ini, kbin, ktbin);
				GraphErrors hBSAy_sub_pi0_sim_tbin_ini = new GraphErrors();
				hBSAy_sub_pi0_sim_tbin.add(hBSAy_sub_pi0_sim_tbin_ini, kbin, ktbin);
				for(int kent = 0; kent < 20; kent++)
				{
					if(bsa_tbin.hasItem(kbin, ktbin, kent) && bsaerr_tbin.hasItem(kbin, ktbin, kent)
							&& phi_rb_poshel_tbin.getItem(kbin, ktbin, kent) + phi_rb_neghel_tbin.getItem(kbin, ktbin, kent) > pentmin)
					{
						hBSAy_tbin.getItem(kbin, ktbin).addPoint(phibc[kent], bsa_tbin.getItem(kbin, ktbin, kent), 0, bsaerr_tbin.getItem(kbin, ktbin, kent));
					}
					if(bsa_tbin_pi0.hasItem(kbin, ktbin, kent) && bsaerr_tbin_pi0.hasItem(kbin, ktbin, kent)
							&& phi_rb_poshel_tbin_pi0.getItem(kbin, ktbin, kent) + phi_rb_neghel_tbin_pi0.getItem(kbin, ktbin, kent) > pentmin_pi0)
					{
						hBSAy_tbin_pi0.getItem(kbin, ktbin).addPoint(phibc[kent], bsa_tbin_pi0.getItem(kbin, ktbin, kent), 0,
																		bsaerr_tbin_pi0.getItem(kbin, ktbin, kent));
					}
					if(f_pi0_sim_tbin.hasItem(kbin, ktbin, kent) && bsa_sub_pi0_sim_tbin.hasItem(kbin, ktbin, kent)
							&& phi_rb_poshel_tbin.getItem(kbin, ktbin, kent) + phi_rb_neghel_tbin.getItem(kbin, ktbin, kent) > pentmin
							&& phi_rb_poshel_tbin_pi0.getItem(kbin, ktbin, kent) + phi_rb_neghel_tbin_pi0.getItem(kbin, ktbin, kent) > pentmin_pi0)
					{
						hf_pi0_sim_tbin.getItem(kbin, ktbin).addPoint(phibc[kent], f_pi0_sim_tbin.getItem(kbin, ktbin, kent), 0,
																		ferr_pi0_sim_tbin.getItem(kbin, ktbin, kent));
						hBSAy_sub_pi0_sim_tbin.getItem(kbin, ktbin).addPoint(phibc[kent], bsa_sub_pi0_sim_tbin.getItem(kbin, ktbin, kent), 0,
																				bsaerr_sub_pi0_sim_tbin.getItem(kbin, ktbin, kent));
					}
				}
			}
		}
		
		JFrame framebsa = new JFrame("BSA");
		framebsa.setSize(1000, 1000);
		EmbeddedCanvas canbsa = new EmbeddedCanvas();
		framebsa.add(canbsa);
		framebsa.setLocationRelativeTo(null);
		framebsa.setVisible(true);
		canbsa.divide(2, 2);
		canbsa.cd(0);
		canbsa.setFont("Arial");
		hBSAy.setTitle("Raw Beam Spin Asymmetry");
		hBSAy.setTitleX("#phi [#degree]");
		hBSAy.setTitleY("Raw BSA");
		hBSAy.setMarkerSize(5);
		hBSAy.setMarkerStyle(1);
		hBSAy.setLineThickness(1);
		canbsa.getPad(0).setTitleFontSize(32);
		canbsa.getPad(0).setAxisTitleFontSize(32);
		canbsa.getPad(0).setAxisLabelFontSize(24);
		canbsa.getPad(0).setAxisRange(0, 360, -0.35, 0.35);
		canbsa.getPad(0).setStatBoxFontSize(18);
		canbsa.draw(hBSAy, "same");
		F1D fBSA= new F1D("fBSA", "[A]*sin(x/57.3)/(1+[B]*cos(x/57.3))", 0, 360);
		fBSA.setParameter(0, 0.2);
		fBSA.setParameter(1, -0.5);
		DataFitter.fit(fBSA, hBSAy, "Q");
		fBSA.setLineColor(2);
		fBSA.setLineWidth(3);
		fBSA.setOptStat(1110);
		canbsa.draw(fBSA, "same");
		canbsa.cd(1);
		canbsa.setFont("Arial");
		hBSAy_pi0.setTitle("DV#pi^0P Beam Spin Asymmetry");
		hBSAy_pi0.setTitleX("#phi [#degree]");
		hBSAy_pi0.setTitleY("DV#pi^0P BSA");
		hBSAy_pi0.setMarkerSize(5);
		hBSAy_pi0.setMarkerStyle(1);
		hBSAy_pi0.setLineThickness(1);
		canbsa.getPad(1).setTitleFontSize(32);
		canbsa.getPad(1).setAxisTitleFontSize(32);
		canbsa.getPad(1).setAxisLabelFontSize(24);
		canbsa.getPad(1).setAxisRange(0, 360, -0.35, 0.35);
		canbsa.getPad(1).setStatBoxFontSize(18);
		canbsa.draw(hBSAy_pi0, "same");
		F1D fBSA_pi0= new F1D("fBSA_pi0", "[A]*sin(x/57.3)", 0, 360);
		fBSA_pi0.setParameter(0, 0.15);
		DataFitter.fit(fBSA_pi0, hBSAy_pi0, "Q");
		fBSA_pi0.setLineColor(2);
		fBSA_pi0.setLineWidth(3);
		fBSA_pi0.setOptStat(1110);
		canbsa.draw(fBSA_pi0, "same");
		canbsa.cd(2);
		canbsa.setFont("Arial");
		hf_pi0_sim.setTitle("#pi_sim f");
		hf_pi0_sim.setTitleX("#phi [#degree]");
		hf_pi0_sim.setTitleY("f");
		hf_pi0_sim.setMarkerSize(5);
		hf_pi0_sim.setMarkerStyle(1);
		hf_pi0_sim.setLineThickness(1);
		canbsa.getPad(2).setTitleFontSize(32);
		canbsa.getPad(2).setAxisTitleFontSize(32);
		canbsa.getPad(2).setAxisLabelFontSize(24);
		canbsa.getPad(2).setAxisRange(0, 360, 0, 1.5);
		canbsa.getPad(2).setStatBoxFontSize(18);
		canbsa.draw(hf_pi0_sim, "same");
		canbsa.cd(3);
		canbsa.setFont("Arial");
		hBSAy_sub_pi0_sim.setTitle("#pi^0_sim Corrected Beam Spin Asymmetry");
		hBSAy_sub_pi0_sim.setTitleX("#phi [#degree]");
		hBSAy_sub_pi0_sim.setTitleY("#pi^0_sim Corrected BSA");
		hBSAy_sub_pi0_sim.setMarkerSize(5);
		hBSAy_sub_pi0_sim.setMarkerStyle(1);
		hBSAy_sub_pi0_sim.setLineThickness(1);
		canbsa.getPad(3).setTitleFontSize(32);
		canbsa.getPad(3).setAxisTitleFontSize(32);
		canbsa.getPad(3).setAxisLabelFontSize(24);
		canbsa.getPad(3).setAxisRange(0, 360, -0.35, 0.35);
		canbsa.getPad(3).setStatBoxFontSize(18);
		canbsa.draw(hBSAy_sub_pi0_sim, "same");
		F1D fBSA_sub_pi0_sim= new F1D("fBSA_sub_pi0_sim", "[A]*sin(x/57.3)/(1+[B]*cos(x/57.3))", 0, 360);
		fBSA_sub_pi0_sim.setParameter(0, 0.2);
		fBSA_sub_pi0_sim.setParameter(1, -0.5);
		DataFitter.fit(fBSA_sub_pi0_sim, hBSAy_sub_pi0_sim, "Q");
		fBSA_sub_pi0_sim.setLineColor(2);
		fBSA_sub_pi0_sim.setLineWidth(3);
		fBSA_sub_pi0_sim.setOptStat(1110);
		canbsa.draw(fBSA_sub_pi0_sim, "same");
		/*
		System.out.println("==========================================================================================");
		System.out.println("BSA:	x	#deltax	y	#delta	y");
		for(int p = 0; p < 20; p++)		
		{
		System.out.println("	" + phibc[p] + "	" + 0 + "	" + bsa.get(p) + "	" + bsaerr.get(p));
		}
		System.out.println("A: " + fBSA.parameter(0).value() + " +/- " + fBSA.parameter(0).error());
		System.out.println("B: " + fBSA.parameter(1).value() + " +/- " + fBSA.parameter(1).error());
		System.out.println("==========================================================================================");
		System.out.println("DV#pi^0P BSA:	x	#deltax	y	#delta	y");
		for(int p = 0; p < 20; p++)		
		{
		System.out.println("	" + phibc[p] + "	" + 0 + "	" + bsa_pi0.get(p) + "	" + bsaerr_pi0.get(p));
		}
		System.out.println("A: " + fBSA_pi0.parameter(0).value() + " +/- " + fBSA_pi0.parameter(0).error());
		System.out.println("==========================================================================================");
		System.out.println("#pi^0_sim f:	x	#deltax	y	#delta	y");
		for(int p = 0; p < 20; p++)		
		{
		System.out.println("	" + phibc[p] + "	" + 0 + "	" + f_pi0_sim.get(p) + "	" + ferr_pi0_sim.get(p));
		}
		*/
		System.out.println("==========================================================================================");
		System.out.println("#pi^0_sim Corrected BSA:	x	y	#delta y");
		for(int p = 0; p < 20; p++)		
		{
		System.out.println("	" + phibc[p] + "	" + bsa_sub_pi0_sim.get(p) + "	" + bsaerr_sub_pi0_sim.get(p));
		}
		System.out.println("A: " + fBSA_sub_pi0_sim.parameter(0).value() + " +/- " + fBSA_sub_pi0_sim.parameter(0).error());
		System.out.println("B: " + fBSA_sub_pi0_sim.parameter(1).value() + " +/- " + fBSA_sub_pi0_sim.parameter(1).error());

		for(int kbin = 0; kbin  < 8; kbin++)
		{
			for(int ktbin = 0; ktbin < 5; ktbin++)
			{
				JFrame framebsa_tbin = new JFrame("BSA: Q^2-x_B Bin " + (kbin+1) + ", -t Bin " + (ktbin+1));
				framebsa_tbin.setSize(1000, 1000);
				EmbeddedCanvas canbsa_tbin = new EmbeddedCanvas();
				framebsa_tbin.add(canbsa_tbin);
				framebsa_tbin.setLocationRelativeTo(null);
				framebsa_tbin.setVisible(true);
				canbsa_tbin.divide(2, 2);
				canbsa_tbin.cd(0);
				canbsa_tbin.setFont("Arial");
				hBSAy_tbin.getItem(kbin, ktbin).setTitle("Raw Beam Spin Asymmetry");
				hBSAy_tbin.getItem(kbin, ktbin).setTitleX("#phi [#degree]");
				hBSAy_tbin.getItem(kbin, ktbin).setTitleY("Raw BSA");
				hBSAy_tbin.getItem(kbin, ktbin).setMarkerSize(5);
				hBSAy_tbin.getItem(kbin, ktbin).setMarkerStyle(1);
				hBSAy_tbin.getItem(kbin, ktbin).setLineThickness(1);
				canbsa_tbin.getPad(0).setTitleFontSize(32);
				canbsa_tbin.getPad(0).setAxisTitleFontSize(32);
				canbsa_tbin.getPad(0).setAxisLabelFontSize(24);
				canbsa_tbin.getPad(0).setAxisRange(0, 360, -0.35, 0.35);
				canbsa_tbin.getPad(0).setStatBoxFontSize(18);
				canbsa_tbin.draw(hBSAy_tbin.getItem(kbin, ktbin), "same");
				F1D ftBSA= new F1D("ftBSA", "[A]*sin(x/57.3)/(1+[B]*cos(x/57.3))", 0, 360);
				ftBSA.setParameter(0, 0.3);
				ftBSA.setParameter(1, -0.5);
				ftBSA.setParLimits(1, -1, 1);
				DataFitter.fit(ftBSA, hBSAy_tbin.getItem(kbin, ktbin), "Q");
				ftBSA.setLineColor(2);
				ftBSA.setLineWidth(3);
				ftBSA.setOptStat(1110);
				canbsa_tbin.draw(ftBSA, "same");
				canbsa_tbin.cd(1);
				canbsa_tbin.setFont("Arial");
				hBSAy_tbin_pi0.getItem(kbin, ktbin).setTitle("DV#pi^0P Beam Spin Asymmetry");
				hBSAy_tbin_pi0.getItem(kbin, ktbin).setTitleX("#phi [#degree]");
				hBSAy_tbin_pi0.getItem(kbin, ktbin).setTitleY("DV#pi^0P BSA");
				hBSAy_tbin_pi0.getItem(kbin, ktbin).setMarkerSize(5);
				hBSAy_tbin_pi0.getItem(kbin, ktbin).setMarkerStyle(1);
				hBSAy_tbin_pi0.getItem(kbin, ktbin).setLineThickness(1);
				canbsa_tbin.getPad(1).setTitleFontSize(32);
				canbsa_tbin.getPad(1).setAxisTitleFontSize(32);
				canbsa_tbin.getPad(1).setAxisLabelFontSize(24);
				canbsa_tbin.getPad(1).setAxisRange(0, 360, -0.35, 0.35);
				canbsa_tbin.getPad(1).setStatBoxFontSize(18);
				canbsa_tbin.draw(hBSAy_tbin_pi0.getItem(kbin, ktbin), "same");
				F1D ftBSA_pi0= new F1D("ftBSA_pi0", "[A]*sin(x/57.3)", 0, 360);
				ftBSA_pi0.setParameter(0, 0.15);
				DataFitter.fit(ftBSA_pi0, hBSAy_tbin_pi0.getItem(kbin, ktbin), "Q");
				ftBSA_pi0.setLineColor(2);
				ftBSA_pi0.setLineWidth(3);
				ftBSA_pi0.setOptStat(110);
				canbsa_tbin.draw(ftBSA_pi0, "same");
				canbsa_tbin.cd(2);
				canbsa_tbin.setFont("Arial");
				hf_pi0_sim_tbin.getItem(kbin, ktbin).setTitle("#pi^0_sim f");
				hf_pi0_sim_tbin.getItem(kbin, ktbin).setTitleX("#phi [#degree]");
				hf_pi0_sim_tbin.getItem(kbin, ktbin).setTitleY("f");
				hf_pi0_sim_tbin.getItem(kbin, ktbin).setMarkerSize(5);
				hf_pi0_sim_tbin.getItem(kbin, ktbin).setMarkerStyle(1);
				hf_pi0_sim_tbin.getItem(kbin, ktbin).setLineThickness(1);
				canbsa_tbin.getPad(2).setTitleFontSize(32);
				canbsa_tbin.getPad(2).setAxisTitleFontSize(32);
				canbsa_tbin.getPad(2).setAxisLabelFontSize(24);
				canbsa_tbin.getPad(2).setAxisRange(0, 360, 0, 1.5);
				canbsa_tbin.getPad(2).setStatBoxFontSize(18);
				canbsa_tbin.draw(hf_pi0_sim_tbin.getItem(kbin, ktbin), "same");
				canbsa_tbin.cd(3);
				canbsa_tbin.setFont("Arial");
				hBSAy_sub_pi0_sim_tbin.getItem(kbin, ktbin).setTitle("#pi^0_sim Corrected Beam Spin Asymmetry");
				hBSAy_sub_pi0_sim_tbin.getItem(kbin, ktbin).setTitleX("#phi [#degree]");
				hBSAy_sub_pi0_sim_tbin.getItem(kbin, ktbin).setTitleY("#pi^0_sim Corrected BSA");
				hBSAy_sub_pi0_sim_tbin.getItem(kbin, ktbin).setMarkerSize(5);
				hBSAy_sub_pi0_sim_tbin.getItem(kbin, ktbin).setMarkerStyle(1);
				hBSAy_sub_pi0_sim_tbin.getItem(kbin, ktbin).setLineThickness(1);
				canbsa_tbin.getPad(3).setTitleFontSize(32);
				canbsa_tbin.getPad(3).setAxisTitleFontSize(32);
				canbsa_tbin.getPad(3).setAxisLabelFontSize(24);
				canbsa_tbin.getPad(3).setAxisRange(0, 360, -0.35, 0.35);
				canbsa_tbin.getPad(3).setStatBoxFontSize(18);
				canbsa_tbin.draw(hBSAy_sub_pi0_sim_tbin.getItem(kbin, ktbin), "same");
				F1D ftBSA_sub_pi0_sim = new F1D("ftBSA_sub_pi0_sim", "[A]*sin(x/57.3)/(1+[B]*cos(x/57.3))", 0, 360);
				ftBSA_sub_pi0_sim.setParameter(0, 0.3);
				ftBSA_sub_pi0_sim.setParameter(1, -0.5);
			//	ftBSA_sub_pi0_sim.setParLimits(1, -1, 1);
				DataFitter.fit(ftBSA_sub_pi0_sim, hBSAy_sub_pi0_sim_tbin.getItem(kbin, ktbin), "Q");
				ftBSA_sub_pi0_sim.setLineColor(2);
				ftBSA_sub_pi0_sim.setLineWidth(3);
				ftBSA_sub_pi0_sim.setOptStat(1110);
				canbsa_tbin.draw(ftBSA_sub_pi0_sim, "same");
				
				if((kbin+1) == 3)
				{
					System.out.println("	");
					System.out.println("BSA: Q^2-x_B Bin " + (kbin+1) + ", -t Bin " + (ktbin+1));
					
					System.out.println("==========================================================================================");
					System.out.println("BSA:	x	N+	N-	y	#delta y");
					for(int p = 0; p < 20; p++)	
					{
						double Np = phi_rb_poshel_tbin.getItem(kbin, ktbin, p);
						double Nm = phi_rb_neghel_tbin.getItem(kbin, ktbin, p);
						double BSA = bsa_tbin.getItem(kbin, ktbin, p);
						double BSAerr = bsaerr_tbin.getItem(kbin, ktbin, p);
						System.out.println("	" + phibc[p] + "	" + Np + "	" + Nm + "	" + BSA + "	" + BSAerr);
					}
					System.out.println("A: " + ftBSA.parameter(0).value() + " +/- " + ftBSA.parameter(0).error());
					System.out.println("B: " + ftBSA.parameter(1).value() + " +/- " + ftBSA.parameter(1).error());
					System.out.println("==========================================================================================");
					System.out.println("DV#pi^0P BSA:	x	N+	N-	y	#delta y");
					for(int p = 0; p < 20; p++)		
					{
						double Np = phi_rb_poshel_tbin_pi0.getItem(kbin, ktbin, p);
						double Nm = phi_rb_neghel_tbin_pi0.getItem(kbin, ktbin, p);
						double BSA = bsa_tbin_pi0.getItem(kbin, ktbin, p);
						double BSAerr = bsaerr_tbin_pi0.getItem(kbin, ktbin, p);
						System.out.println("	" + phibc[p] + "	" + Np + "	" + Nm + "	" + BSA + "	" + BSAerr);
					}
					System.out.println("A: " + ftBSA_pi0.parameter(0).value() + " +/- " + ftBSA_pi0.parameter(0).error());
					System.out.println("==========================================================================================");
					System.out.println("#pi^0 f:	x	Ng	Ngg	y	#delta y");
					for(int p = 0; p < 20; p++)		
					{
						double Ng = phi_rb_1gamma_pi0_sim_tbin.getItem(kbin, ktbin, p);
						double Ngg = phi_rb_2gamma_pi0_sim_tbin.getItem(kbin, ktbin, p);
						double f = f_pi0_sim_tbin.getItem(kbin, ktbin, p);
						double ferr = ferr_pi0_sim_tbin.getItem(kbin, ktbin, p);
						System.out.println("	" + phibc[p] + "	" + Ng + "	" + Ngg
										+ "	" + f + "	" + ferr);
					}
					
					System.out.println("==========================================================================================");
					System.out.println("#pi0_sim Corrected BSA:	x	#deltax	y	#delta	y");
					for(int p = 0; p < 20; p++)		
					{
					System.out.println("	" + phibc[p] + "	" + 0 + "	" + bsa_sub_pi0_sim_tbin.getItem(kbin, ktbin, p)
										+ "	" + bsaerr_sub_pi0_sim_tbin.getItem(kbin, ktbin, p));
					}
					System.out.println("A: " + ftBSA_sub_pi0_sim.parameter(0).value() + " +/- " + ftBSA_sub_pi0_sim.parameter(0).error());
					System.out.println("B: " + ftBSA_sub_pi0_sim.parameter(1).value() + " +/- " + ftBSA_sub_pi0_sim.parameter(1).error());
				}
				/*
				System.out.println("	");
				System.out.println("BSA: Q^2-x_B Bin " + (kbin+1) + ", -t Bin " + (ktbin+1));
				
				System.out.println("==========================================================================================");
				System.out.println("BSA:	x	#deltax	y	#delta	y");
				for(int p = 0; p < 20; p++)		
				{
				System.out.println("	" + phibc[p] + "	" + 0 + "	" + bsa_tbin.getItem(kbin, ktbin, p)
									+ "	" + bsaerr_tbin.getItem(kbin, ktbin, p));
				}
				System.out.println("A: " + ftBSA.parameter(0).value() + " +/- " + ftBSA.parameter(0).error());
				System.out.println("B: " + ftBSA.parameter(1).value() + " +/- " + ftBSA.parameter(1).error());
				System.out.println("==========================================================================================");
				System.out.println("DV#pi^0P BSA:	x	#deltax	y	#delta	y");
				for(int p = 0; p < 20; p++)		
				{
				System.out.println("	" + phibc[p] + "	" + 0 + "	" + bsa_tbin_pi0.getItem(kbin, ktbin, p)
									+ "	" + bsaerr_tbin_pi0.getItem(kbin, ktbin, p));
				}
				System.out.println("A: " + ftBSA_pi0.parameter(0).value() + " +/- " + ftBSA_pi0.parameter(0).error());
				System.out.println("==========================================================================================");
				System.out.println("#pi^0 f:	x	#deltax	y	#delta	y");
				for(int p = 0; p < 20; p++)		
				{
				System.out.println("	" + phibc[p] + "	" + 0 + "	" + f_pi0_sim_tbin.getItem(kbin, ktbin, p)
									+ "	" + ferr_pi0_sim_tbin.getItem(kbin, ktbin, p));
				}
				
				System.out.println("==========================================================================================");
				System.out.println("#pi0_sim Corrected BSA:	x	#deltax	y	#delta	y");
				for(int p = 0; p < 20; p++)		
				{
				System.out.println("	" + phibc[p] + "	" + 0 + "	" + bsa_sub_pi0_sim_tbin.getItem(kbin, ktbin, p)
									+ "	" + bsaerr_sub_pi0_sim_tbin.getItem(kbin, ktbin, p));
				}
				*/
				System.out.println("A: " + ftBSA_sub_pi0_sim.parameter(0).value() + " +/- " + ftBSA_sub_pi0_sim.parameter(0).error());
				System.out.println("B: " + ftBSA_sub_pi0_sim.parameter(1).value() + " +/- " + ftBSA_sub_pi0_sim.parameter(1).error());
				halpha_vs_t.getItem(kbin).addPoint(tba.getItem(kbin)[ktbin], ftBSA.parameter(0).value(), 0, ftBSA.parameter(0).error());
				halpha_vs_t_pi0.getItem(kbin).addPoint(tba.getItem(kbin)[ktbin], ftBSA_pi0.parameter(0).value(), 0, ftBSA_pi0.parameter(0).error());
				halpha_corr_pi0_sim_vs_t.getItem(kbin).addPoint(tba.getItem(kbin)[ktbin], ftBSA_sub_pi0_sim.parameter(0).value(), 0, ftBSA_sub_pi0_sim.parameter(0).error());
			}
		}
		
		JFrame frame_alphat = new JFrame("#alpha vs. -t");
		frame_alphat.setSize(1000, 1000);
		EmbeddedCanvas can_alphat = new EmbeddedCanvas();
		frame_alphat.add(can_alphat);
		frame_alphat.setLocationRelativeTo(null);
		frame_alphat.setVisible(true);
		can_alphat.divide(2, 2);
		can_alphat.cd(0);
		can_alphat.setFont("Arial");
		halpha_vs_t.getItem(0).setTitle("Raw BSA #alpha vs. -t");
		halpha_vs_t.getItem(0).setTitleX("-t [GeV^2]");
		halpha_vs_t.getItem(0).setTitleY("Pre-Correction #alpha");
		halpha_vs_t.getItem(0).setMarkerSize(5);
		halpha_vs_t.getItem(0).setMarkerStyle(1);
		halpha_vs_t.getItem(0).setLineThickness(1);
		can_alphat.getPad(0).setTitleFontSize(32);
		can_alphat.getPad(0).setAxisTitleFontSize(32);
		can_alphat.getPad(0).setAxisLabelFontSize(24);
		can_alphat.getPad(0).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat.getPad(0).setStatBoxFontSize(18);
		can_alphat.draw(halpha_vs_t.getItem(0), "same");
		can_alphat.cd(1);
		can_alphat.setFont("Arial");
		halpha_vs_t.getItem(1).setTitle("Raw BSA #alpha vs. -t");
		halpha_vs_t.getItem(1).setTitleX("-t [GeV^2]");
		halpha_vs_t.getItem(1).setTitleY("Pre-Correction #alpha");
		halpha_vs_t.getItem(1).setMarkerSize(5);
		halpha_vs_t.getItem(1).setMarkerStyle(1);
		halpha_vs_t.getItem(1).setLineThickness(1);
		can_alphat.getPad(1).setTitleFontSize(32);
		can_alphat.getPad(1).setAxisTitleFontSize(32);
		can_alphat.getPad(1).setAxisLabelFontSize(24);
		can_alphat.getPad(1).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat.getPad(1).setStatBoxFontSize(18);
		can_alphat.draw(halpha_vs_t.getItem(1), "same");
		halpha_vs_t.getItem(2).setMarkerSize(5);
		halpha_vs_t.getItem(2).setLineThickness(1);
		halpha_vs_t.getItem(2).setMarkerStyle(0);
		halpha_vs_t.getItem(2).setMarkerColor(2);
		halpha_vs_t.getItem(2).setLineColor(2);
		can_alphat.draw(halpha_vs_t.getItem(2), "same");
		can_alphat.cd(2);
		can_alphat.setFont("Arial");
		halpha_vs_t.getItem(3).setTitle("Raw BSA #alpha vs. -t");
		halpha_vs_t.getItem(3).setTitleX("-t [GeV^2]");
		halpha_vs_t.getItem(3).setTitleY("Pre-Correction #alpha");
		halpha_vs_t.getItem(3).setMarkerSize(5);
		halpha_vs_t.getItem(3).setMarkerStyle(1);
		halpha_vs_t.getItem(3).setLineThickness(1);
		can_alphat.getPad(2).setTitleFontSize(32);
		can_alphat.getPad(2).setAxisTitleFontSize(32);
		can_alphat.getPad(2).setAxisLabelFontSize(24);
		can_alphat.getPad(2).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat.getPad(2).setStatBoxFontSize(18);
		can_alphat.draw(halpha_vs_t.getItem(3), "same");
		halpha_vs_t.getItem(4).setMarkerSize(5);
		halpha_vs_t.getItem(4).setLineThickness(1);
		halpha_vs_t.getItem(4).setMarkerStyle(0);
		halpha_vs_t.getItem(4).setMarkerColor(2);
		halpha_vs_t.getItem(4).setLineColor(2);
		can_alphat.draw(halpha_vs_t.getItem(4), "same");
		halpha_vs_t.getItem(5).setMarkerSize(5);
		halpha_vs_t.getItem(5).setLineThickness(1);
		halpha_vs_t.getItem(5).setMarkerStyle(2);
		halpha_vs_t.getItem(5).setMarkerColor(4);
		halpha_vs_t.getItem(5).setLineColor(4);
		can_alphat.draw(halpha_vs_t.getItem(5), "same");
		can_alphat.cd(3);
		can_alphat.setFont("Arial");
		halpha_vs_t.getItem(6).setTitle("Raw BSA #alpha vs. -t'");
		halpha_vs_t.getItem(6).setTitleX("-t [GeV^2]");
		halpha_vs_t.getItem(6).setTitleY("Pre-Correction #alpha");
		halpha_vs_t.getItem(6).setMarkerSize(5);
		halpha_vs_t.getItem(6).setMarkerStyle(1);
		halpha_vs_t.getItem(6).setLineThickness(1);
		can_alphat.getPad(3).setTitleFontSize(32);
		can_alphat.getPad(3).setAxisTitleFontSize(32);
		can_alphat.getPad(3).setAxisLabelFontSize(24);
		can_alphat.getPad(3).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat.getPad(3).setStatBoxFontSize(18);
		can_alphat.draw(halpha_vs_t.getItem(6), "same");
		halpha_vs_t.getItem(7).setMarkerSize(5);
		halpha_vs_t.getItem(7).setLineThickness(1);
		halpha_vs_t.getItem(7).setMarkerStyle(0);
		halpha_vs_t.getItem(7).setMarkerColor(2);
		halpha_vs_t.getItem(7).setLineColor(2);
		can_alphat.draw(halpha_vs_t.getItem(7), "same");
		
		JFrame frame_alphat_pi0 = new JFrame("#pi^0 #alpha vs. -t");
		frame_alphat_pi0.setSize(1000, 1000);
		EmbeddedCanvas can_alphat_pi0 = new EmbeddedCanvas();
		frame_alphat_pi0.add(can_alphat_pi0);
		frame_alphat_pi0.setLocationRelativeTo(null);
		frame_alphat_pi0.setVisible(true);
		can_alphat_pi0.divide(2, 2);
		can_alphat_pi0.cd(0);
		can_alphat_pi0.setFont("Arial");
		halpha_vs_t_pi0.getItem(0).setTitle("DV#pi^0P BSA #alpha vs. -t");
		halpha_vs_t_pi0.getItem(0).setTitleX("-t [GeV^2]");
		halpha_vs_t_pi0.getItem(0).setTitleY("DV#pi^0P #alpha");
		halpha_vs_t_pi0.getItem(0).setMarkerSize(5);
		halpha_vs_t_pi0.getItem(0).setMarkerStyle(1);
		halpha_vs_t_pi0.getItem(0).setLineThickness(1);
		can_alphat_pi0.getPad(0).setTitleFontSize(32);
		can_alphat_pi0.getPad(0).setAxisTitleFontSize(32);
		can_alphat_pi0.getPad(0).setAxisLabelFontSize(24);
		can_alphat_pi0.getPad(0).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_pi0.getPad(0).setStatBoxFontSize(18);
		can_alphat_pi0.draw(halpha_vs_t_pi0.getItem(0), "same");
		can_alphat_pi0.cd(1);
		can_alphat_pi0.setFont("Arial");
		halpha_vs_t_pi0.getItem(1).setTitle("DV#pi^0P BSA #alpha vs. -t");
		halpha_vs_t_pi0.getItem(1).setTitleX("-t [GeV^2]");
		halpha_vs_t_pi0.getItem(1).setTitleY("DV#pi^0P #alpha");
		halpha_vs_t_pi0.getItem(1).setMarkerSize(5);
		halpha_vs_t_pi0.getItem(1).setMarkerStyle(1);
		halpha_vs_t_pi0.getItem(1).setLineThickness(1);
		can_alphat_pi0.getPad(1).setTitleFontSize(32);
		can_alphat_pi0.getPad(1).setAxisTitleFontSize(32);
		can_alphat_pi0.getPad(1).setAxisLabelFontSize(24);
		can_alphat_pi0.getPad(1).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_pi0.getPad(1).setStatBoxFontSize(18);
		can_alphat_pi0.draw(halpha_vs_t_pi0.getItem(1), "same");
		halpha_vs_t_pi0.getItem(2).setMarkerSize(5);
		halpha_vs_t_pi0.getItem(2).setLineThickness(1);
		halpha_vs_t_pi0.getItem(2).setMarkerStyle(0);
		halpha_vs_t_pi0.getItem(2).setMarkerColor(2);
		halpha_vs_t_pi0.getItem(2).setLineColor(2);
		can_alphat_pi0.draw(halpha_vs_t_pi0.getItem(2), "same");
		can_alphat_pi0.cd(2);
		can_alphat_pi0.setFont("Arial");
		halpha_vs_t_pi0.getItem(3).setTitle("DV#pi^0P BSA #alpha vs. -t");
		halpha_vs_t_pi0.getItem(3).setTitleX("-t [GeV^2]");
		halpha_vs_t_pi0.getItem(3).setTitleY("DV#pi^0P #alpha");
		halpha_vs_t_pi0.getItem(3).setMarkerSize(5);
		halpha_vs_t_pi0.getItem(3).setMarkerStyle(1);
		halpha_vs_t_pi0.getItem(3).setLineThickness(1);
		can_alphat_pi0.getPad(2).setTitleFontSize(32);
		can_alphat_pi0.getPad(2).setAxisTitleFontSize(32);
		can_alphat_pi0.getPad(2).setAxisLabelFontSize(24);
		can_alphat_pi0.getPad(2).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_pi0.getPad(2).setStatBoxFontSize(18);
		can_alphat_pi0.draw(halpha_vs_t_pi0.getItem(3), "same");
		halpha_vs_t_pi0.getItem(4).setMarkerSize(5);
		halpha_vs_t_pi0.getItem(4).setLineThickness(1);
		halpha_vs_t_pi0.getItem(4).setMarkerStyle(0);
		halpha_vs_t_pi0.getItem(4).setMarkerColor(2);
		halpha_vs_t_pi0.getItem(4).setLineColor(2);
		can_alphat_pi0.draw(halpha_vs_t_pi0.getItem(4), "same");
		halpha_vs_t_pi0.getItem(5).setMarkerSize(5);
		halpha_vs_t_pi0.getItem(5).setLineThickness(1);
		halpha_vs_t_pi0.getItem(5).setMarkerStyle(2);
		halpha_vs_t_pi0.getItem(5).setMarkerColor(4);
		halpha_vs_t_pi0.getItem(5).setLineColor(4);
		can_alphat_pi0.draw(halpha_vs_t_pi0.getItem(5), "same");
		can_alphat_pi0.cd(3);
		can_alphat_pi0.setFont("Arial");
		halpha_vs_t_pi0.getItem(6).setTitle("DV#pi^0P BSA #alpha vs. -t'");
		halpha_vs_t_pi0.getItem(6).setTitleX("-t [GeV^2]");
		halpha_vs_t_pi0.getItem(6).setTitleY("DV#pi^0P #alpha");
		halpha_vs_t_pi0.getItem(6).setMarkerSize(5);
		halpha_vs_t_pi0.getItem(6).setMarkerStyle(1);
		halpha_vs_t_pi0.getItem(6).setLineThickness(1);
		can_alphat_pi0.getPad(3).setTitleFontSize(32);
		can_alphat_pi0.getPad(3).setAxisTitleFontSize(32);
		can_alphat_pi0.getPad(3).setAxisLabelFontSize(24);
		can_alphat_pi0.getPad(3).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_pi0.getPad(3).setStatBoxFontSize(18);
		can_alphat_pi0.draw(halpha_vs_t_pi0.getItem(6), "same");
		halpha_vs_t_pi0.getItem(7).setMarkerSize(5);
		halpha_vs_t_pi0.getItem(7).setLineThickness(1);
		halpha_vs_t_pi0.getItem(7).setMarkerStyle(0);
		halpha_vs_t_pi0.getItem(7).setMarkerColor(2);
		halpha_vs_t_pi0.getItem(7).setLineColor(2);
		can_alphat_pi0.draw(halpha_vs_t_pi0.getItem(7), "same");
		
		JFrame frame_alphat_corr_pi0_sim = new JFrame("#pi^0_sim Corrected #alpha vs. -t");
		frame_alphat_corr_pi0_sim.setSize(1000, 1000);
		EmbeddedCanvas can_alphat_corr_pi0_sim = new EmbeddedCanvas();
		frame_alphat_corr_pi0_sim.add(can_alphat_corr_pi0_sim);
		frame_alphat_corr_pi0_sim.setLocationRelativeTo(null);
		frame_alphat_corr_pi0_sim.setVisible(true);
		can_alphat_corr_pi0_sim.divide(2, 2);
		can_alphat_corr_pi0_sim.cd(0);
		can_alphat_corr_pi0_sim.setFont("Arial");
		halpha_corr_pi0_sim_vs_t.getItem(0).setTitle("#pi^0_sim Corrected BSA #alpha vs. -t");
		halpha_corr_pi0_sim_vs_t.getItem(0).setTitleX("-t [GeV^2]");
		halpha_corr_pi0_sim_vs_t.getItem(0).setTitleY("#pi^0_sim Corrected #alpha");
		halpha_corr_pi0_sim_vs_t.getItem(0).setMarkerSize(5);
		halpha_corr_pi0_sim_vs_t.getItem(0).setMarkerStyle(1);
		halpha_corr_pi0_sim_vs_t.getItem(0).setLineThickness(1);
		can_alphat_corr_pi0_sim.getPad(0).setTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(0).setAxisTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(0).setAxisLabelFontSize(24);
		can_alphat_corr_pi0_sim.getPad(0).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_corr_pi0_sim.getPad(0).setStatBoxFontSize(18);
		can_alphat_corr_pi0_sim.draw(halpha_corr_pi0_sim_vs_t.getItem(0), "same");
		can_alphat_corr_pi0_sim.cd(1);
		can_alphat_corr_pi0_sim.setFont("Arial");
		halpha_corr_pi0_sim_vs_t.getItem(1).setTitle("#pi^0_sim Corrected BSA #alpha vs. -t");
		halpha_corr_pi0_sim_vs_t.getItem(1).setTitleX("-t [GeV^2]");
		halpha_corr_pi0_sim_vs_t.getItem(1).setTitleY("#pi^0_sim Corrected #alpha");
		halpha_corr_pi0_sim_vs_t.getItem(1).setMarkerSize(5);
		halpha_corr_pi0_sim_vs_t.getItem(1).setMarkerStyle(1);
		halpha_corr_pi0_sim_vs_t.getItem(1).setLineThickness(1);
		can_alphat_corr_pi0_sim.getPad(1).setTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(1).setAxisTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(1).setAxisLabelFontSize(24);
		can_alphat_corr_pi0_sim.getPad(1).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_corr_pi0_sim.getPad(1).setStatBoxFontSize(18);
		can_alphat_corr_pi0_sim.draw(halpha_corr_pi0_sim_vs_t.getItem(1), "same");
		halpha_corr_pi0_sim_vs_t.getItem(2).setMarkerSize(5);
		halpha_corr_pi0_sim_vs_t.getItem(2).setLineThickness(1);
		halpha_corr_pi0_sim_vs_t.getItem(2).setMarkerStyle(0);
		halpha_corr_pi0_sim_vs_t.getItem(2).setMarkerColor(2);
		halpha_corr_pi0_sim_vs_t.getItem(2).setLineColor(2);
		can_alphat_corr_pi0_sim.draw(halpha_corr_pi0_sim_vs_t.getItem(2), "same");
		can_alphat_corr_pi0_sim.cd(2);
		can_alphat_corr_pi0_sim.setFont("Arial");
		halpha_corr_pi0_sim_vs_t.getItem(3).setTitle("#pi^0_sim Corrected BSA #alpha vs. -t");
		halpha_corr_pi0_sim_vs_t.getItem(3).setTitleX("-t [GeV^2]");
		halpha_corr_pi0_sim_vs_t.getItem(3).setTitleY("#pi^0_sim Corrected #alpha");
		halpha_corr_pi0_sim_vs_t.getItem(3).setMarkerSize(5);
		halpha_corr_pi0_sim_vs_t.getItem(3).setMarkerStyle(1);
		halpha_corr_pi0_sim_vs_t.getItem(3).setLineThickness(1);
		can_alphat_corr_pi0_sim.getPad(2).setTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(2).setAxisTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(2).setAxisLabelFontSize(24);
		can_alphat_corr_pi0_sim.getPad(2).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_corr_pi0_sim.getPad(2).setStatBoxFontSize(18);
		can_alphat_corr_pi0_sim.draw(halpha_corr_pi0_sim_vs_t.getItem(3), "same");
		halpha_corr_pi0_sim_vs_t.getItem(4).setMarkerSize(5);
		halpha_corr_pi0_sim_vs_t.getItem(4).setLineThickness(1);
		halpha_corr_pi0_sim_vs_t.getItem(4).setMarkerStyle(0);
		halpha_corr_pi0_sim_vs_t.getItem(4).setMarkerColor(2);
		halpha_corr_pi0_sim_vs_t.getItem(4).setLineColor(2);
		can_alphat_corr_pi0_sim.draw(halpha_corr_pi0_sim_vs_t.getItem(4), "same");
		halpha_corr_pi0_sim_vs_t.getItem(5).setMarkerSize(5);
		halpha_corr_pi0_sim_vs_t.getItem(5).setLineThickness(1);
		halpha_corr_pi0_sim_vs_t.getItem(5).setMarkerStyle(2);
		halpha_corr_pi0_sim_vs_t.getItem(5).setMarkerColor(4);
		halpha_corr_pi0_sim_vs_t.getItem(5).setLineColor(4);
		can_alphat_corr_pi0_sim.draw(halpha_corr_pi0_sim_vs_t.getItem(5), "same");
		can_alphat_corr_pi0_sim.cd(3);
		can_alphat_corr_pi0_sim.setFont("Arial");
		halpha_corr_pi0_sim_vs_t.getItem(6).setTitle("#pi^0_sim Corrected BSA #alpha vs. -t'");
		halpha_corr_pi0_sim_vs_t.getItem(6).setTitleX("-t [GeV^2]");
		halpha_corr_pi0_sim_vs_t.getItem(6).setTitleY("#pi^0_sim Corrected #alpha");
		halpha_corr_pi0_sim_vs_t.getItem(6).setMarkerSize(5);
		halpha_corr_pi0_sim_vs_t.getItem(6).setMarkerStyle(1);
		halpha_corr_pi0_sim_vs_t.getItem(6).setLineThickness(1);
		can_alphat_corr_pi0_sim.getPad(3).setTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(3).setAxisTitleFontSize(32);
		can_alphat_corr_pi0_sim.getPad(3).setAxisLabelFontSize(24);
		can_alphat_corr_pi0_sim.getPad(3).setAxisRange(0, 2, -0.1, 0.45);
		can_alphat_corr_pi0_sim.getPad(3).setStatBoxFontSize(18);
		can_alphat_corr_pi0_sim.draw(halpha_corr_pi0_sim_vs_t.getItem(6), "same");
		halpha_corr_pi0_sim_vs_t.getItem(7).setMarkerSize(5);
		halpha_corr_pi0_sim_vs_t.getItem(7).setLineThickness(1);
		halpha_corr_pi0_sim_vs_t.getItem(7).setMarkerStyle(0);
		halpha_corr_pi0_sim_vs_t.getItem(7).setMarkerColor(2);
		halpha_corr_pi0_sim_vs_t.getItem(7).setLineColor(2);
		can_alphat_corr_pi0_sim.draw(halpha_corr_pi0_sim_vs_t.getItem(7), "same");
	}
}