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

public class Q2_xB_t_binner {
	
	static HipoDataSource reader = new HipoDataSource();
	
	static H2F hQ2_vs_x_B_rb = new H2F("Q2_vs_x_B_rb", "Q2_vs_x_B_rb", 100, 0, 0.7, 100, 0, 7);
	static H2F hnegt_vs_phi_trento = new H2F("-t_vs_phi_trento", "-t_vs_phi_trento", 100, 0, 360, 100, 0, 3);
	
	static H1F hx_B_rb = new H1F("x_B_rb", "x_B_rb", 1000, 0, 0.8);
	static H1F hnegt = new H1F("-t", "-t", 1000, 0, 3);
	
	static IndexedList<H1F> histGroups_x_B_bin = new IndexedList<H1F>(1);
	public static void x_B_bin_histos() {
		for(int ibin = 0; ibin < 4; ibin++) {
			H1F x_B_bin = new H1F("x_B_bin", "x_B_bin", 1000, 0, 0.8);
			histGroups_x_B_bin.add(x_B_bin, ibin);
		}
	}
	
	static IndexedList<H1F> histGroups_x_B_bin_Q2 = new IndexedList<H1F>(1);
	public static void x_B_bin_Q2_histos() {
		for(int ibin = 0; ibin < 4; ibin++) {
			H1F x_B_bin_Q2 = new H1F("x_B_bin_Q2", "x_B_bin_Q2", 1000, 0, 8);
			histGroups_x_B_bin_Q2.add(x_B_bin_Q2, ibin);
		}
	}
	
	static IndexedList<H1F> histGroups_Q2_x_B_bin = new IndexedList<H1F>(1);
	public static void Q2_x_B_bin_histos() {
		for(int ibin = 0; ibin < 8; ibin++) {
			H1F Q2_x_B_bin = new H1F("Q2_x_B_bin", "Q2_x_B_bin", 1000, 0, 8);
			histGroups_Q2_x_B_bin.add(Q2_x_B_bin, ibin);
		}
	}
	
	static IndexedList<H2F> histGroups_negt_vs_phi_trento_bin = new IndexedList<H2F>(1);
	public static void negt_vs_phi_trento_bin_histos() {
		for(int ibin = 0; ibin < 8; ibin++) {
			H2F negt_vs_phi_trento_bin = new H2F("negt_vs_phi_trento_bin", "negt_vs_phi_trento_bin", 100, 0, 360, 100, 0, 3);
			histGroups_negt_vs_phi_trento_bin.add(negt_vs_phi_trento_bin, ibin);
		}
	}
	
	static IndexedList<H1F> histGroups_negt_bin = new IndexedList<H1F>(1);
	public static void negt_bin_histos() {
		for(int ibin = 0; ibin < 5; ibin++) {
			H1F negt_bin = new H1F("negt_bin", "negt_bin", 1000, 0, 3);
			histGroups_negt_bin.add(negt_bin, ibin);
		}
	}
	
	static IndexedList<H1F> histGroups_Q2_x_B_bin_negt = new IndexedList<H1F>(1);
	public static void Q2_x_B_bin_negt_histos() {
		for(int ibin = 0; ibin < 8; ibin++) {
			H1F Q2_x_B_bin_negt = new H1F("Q2_x_B_bin_negt", "Q2_x_B_bin_negt", 1000, 0, 3);
			histGroups_Q2_x_B_bin_negt.add(Q2_x_B_bin_negt, ibin);
		}
	}
	
	static IndexedList<H1F> histGroups_Q2_x_B_negt_bin = new IndexedList<H1F>(2);
	public static void Q2_x_B_negt_bin_histos() {
		for(int ibin = 0; ibin < 8; ibin++) {
			for(int itbin = 0; itbin < 5; itbin++)
			{
				H1F Q2_x_B_negt_bin = new H1F("Q2_x_B_negt_bin", "Q2_x_B_negt_bin", 1000, 0, 3);
				histGroups_Q2_x_B_negt_bin.add(Q2_x_B_negt_bin, ibin, itbin);
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
			int X0count = 0;
			int e_pindex = -1;
			boolean fid_ecal_e = false;
			int proton_pindex = -1;
			int gamma1_pindex = -1;
			int gamma2_pindex = -1;
			int X0_pindex = -1;
			float beta_gamma = 0;
			float beta_X0 = 0;
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0.000511);
			LorentzVector p_proton = new LorentzVector(0, 0, 0, 0.938272);
			LorentzVector p_gamma1 = new LorentzVector(0, 0, 0, 0);
			LorentzVector p_gamma2 = new LorentzVector(0, 0, 0, 0);
			LorentzVector p_X0 = new LorentzVector(0, 0, 0, 0);
// Marker 1
		//	double E = 6.535;
			double E = 7.546;
			double E_gamma_cut = 0.3;
			float calE1 = 0;
			float calE2 = 0;
			float pcalv_e = 0;
			float pcalw_e = 0;
			byte det_gamma1 = 0;
			float pcalv_gamma1 = 0;
			float pcalw_gamma1 = 0;
			float ftcalrad_gamma1 = 0;
			byte det_gamma2 = 0;
			float pcalv_gamma2 = 0;
			float pcalw_gamma2 = 0;
			float ftcalrad_gamma2 = 0;
			byte calsector = 0;
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
					beta_gamma = beta;
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
			
// Marker 2
			
		// DVCS
			if(ecount == 1 && protoncount == 1 && gammacount == 1
					&& v_e.pz() > -12 && v_e.pz() < 7 && Math.abs(v_e.pz()-v_proton.pz()) < 2.5 + (2.5/p_proton.p())
					&& fid_ecal_e == true
					&& (fid_ecal_gamma1 == true || fid_ftcal_gamma1 == true)
					&& p_gamma1.vect().theta(p_e.vect()) > 5
					&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30))
					&& 57.3*p_proton.theta() < 75)
			{
			
				
		// 7.5 GeV DVCS
				if(((det_gamma1 == 10
						&& p_gamma1.vect().theta(X_gamma1.vect()) < 3
						&& Math.sqrt((X.px()*X.px())+(X.py()*X.py())) < 0.25
						&& X.e() < 0.75
						&& X_proton.mass() > 0.5 && X_proton.mass() < 1.6)
						|| (det_gamma1 == 7
						&& p_gamma1.vect().theta(X_gamma1.vect()) < 3.5
						&& Math.sqrt((X.px()*X.px())+(X.py()*X.py())) < 0.4
						&& X.e() < 1.5
						&& X_proton.mass() > 0.2 && X_proton.mass() < 1.9))
						&& Q2 > 1 && W > 2
						&& -t_cal.mass2() > 0.938272*0.938272*x_B*x_B/(1-x_B+x_B*0.938272*0.938272/Q2*Q2))
				
				/*
		// 6.5 GeV DVCS
				if((det_gamma1 == 7
						&& p_gamma1.vect().theta(X_gamma1.vect()) < 4.5
						&& Math.sqrt((X.px()*X.px())+(X.py()*X.py())) < 0.3
						&& X.e() < 1.25
						&& X_proton.mass() > 0.4 && X_proton.mass() < 1.8)
						&& Q2 > 1 && W > 2
						&& -t_cal.mass2() > 0.938272*0.938272*x_B*x_B/(1-x_B+x_B*0.938272*0.938272/Q2*Q2))
				*/
			/*
	// pi^0	
			if(ecount == 1 && protoncount == 1 && gammacount == 2
				&& p_gamma1.vect().theta(p_gamma2.vect()) > 2
				&& v_e.pz() > -12 && v_e.pz() < 7 && Math.abs(v_e.pz()-v_proton.pz()) < 2.5 + (2.5/p_proton.p())
				&& pcalv_e > 12 && pcalw_e > 12
				&& (( pcalv_gamma1 > 12 && pcalw_gamma1 > 12) || (ftcalrad_gamma1 > 8 && ftcalrad_gamma1 < 15))
				&& (( pcalv_gamma2 > 12 && pcalw_gamma2 > 12) || (ftcalrad_gamma2 > 8 && ftcalrad_gamma2 < 15))
				&& p_gamma1.vect().theta(p_e.vect()) > 5
				&& p_gamma2.vect().theta(p_e.vect()) > 5
				&& p_gamma1.vect().theta(p_gamma2.vect()) > 2
				&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30))
				&& 57.3*p_proton.theta() < 75)
			{
				X_proton.sub(p_gamma2);
				X.sub(p_gamma2);
				t_cal.sub(p_gamma2);
				p_gamma1.add(p_gamma2);
			*/
				/*
	// 6.5 GeV pi^0
				if((det_gamma1 == 7 && det_gamma2 == 7
						&& p_gamma1.vect().theta(X_gamma1.vect()) < 4.5
						&& Math.sqrt((X.px()*X.px())+(X.py()*X.py())) < 0.3
						&& X.e() < 1.25 && X_proton.mass() > 0.4 && X_proton.mass() < 1.8)
						&& Q2 > 1 && W > 2
						&& p_gamma1.mass() > 0.103 && p_gamma1.mass() < 0.163
						&& p_gamma1.vect().theta(p_gamma2.vect()) > 2.5)
				*/
				{	
					hQ2_vs_x_B_rb.fill(x_B, Q2);
					hnegt_vs_phi_trento.fill(phi, -t_cal.mass2());
					hx_B_rb.fill(x_B);
					
					if(x_B < 0.12) histGroups_x_B_bin_Q2.getItem(0).fill(Q2);
					else if(x_B >= 0.12 && x_B < 0.16) histGroups_x_B_bin_Q2.getItem(1).fill(Q2);
					else if(x_B >= 0.16 && x_B < 0.24) histGroups_x_B_bin_Q2.getItem(2).fill(Q2);
					else if(x_B >= 0.24) histGroups_x_B_bin_Q2.getItem(3).fill(Q2);
					
					int Q2xBbin = -1;
					if(x_B < 0.12) Q2xBbin = 0;
					else if(x_B >= 0.12 && x_B < 0.16 && Q2 < 1.3) Q2xBbin = 1;
					else if(x_B >= 0.12 && x_B < 0.16 && Q2 >= 1.3) Q2xBbin = 2;
					else if(x_B >= 0.16 && x_B < 0.24 && Q2 < 1.3) Q2xBbin = 3;
					else if(x_B >= 0.16 && x_B < 0.24 && Q2 >= 1.3 && Q2 < 1.7) Q2xBbin = 4;
					else if(x_B >= 0.16 && x_B < 0.24 && Q2 >= 1.7) Q2xBbin = 5;
					else if(x_B >= 0.24 && Q2 < 2.0) Q2xBbin = 6;
					else if(x_B >= 0.24 && Q2 >= 2.0) Q2xBbin = 7;
					
					histGroups_negt_vs_phi_trento_bin.getItem(Q2xBbin).fill(phi, -t_cal.mass2());
					histGroups_Q2_x_B_bin_negt.getItem(Q2xBbin).fill(-t_cal.mass2());
					
// Marker 3					
					/* if(x_B >= 0.12 && x_B < 0.16) */ hnegt.fill(-t_cal.mass2());
				}
			}
		}
	}
	
	public static void main(String[] args) {
		
		x_B_bin_histos();
		x_B_bin_Q2_histos();
		Q2_x_B_bin_histos();
		negt_vs_phi_trento_bin_histos();
		negt_bin_histos();
		Q2_x_B_bin_negt_histos();
		Q2_x_B_negt_bin_histos();

// Marker 4
		reader.open("C:/Users/joshtanj/Documents/download/skim_epg_bank_merged_7546MeV_skim16.hipo");
		
		int eventCounter = 1;
		while(reader.hasEvent() && eventCounter < 10000000)
		{
			processEvent(reader.getNextEvent(), eventCounter);
			if(eventCounter%50000 == 0) System.out.println("Event: " + eventCounter);
			eventCounter++;
		}
		int Nevent = eventCounter-1;
		System.out.println("Number of events: " + Nevent);
		
		JFrame framekin = new JFrame("Kinematic Variables");
		framekin.setSize(1000, 500);
		EmbeddedCanvas cankin = new EmbeddedCanvas();
		framekin.add(cankin);
		framekin.setLocationRelativeTo(null);
		framekin.setVisible(true);
		cankin.divide(2, 1);
		cankin.cd(0);
		cankin.setFont("Arial");
		hQ2_vs_x_B_rb.setTitle("Q^2 vs. x_B");
		hQ2_vs_x_B_rb.setTitleX("x_B");
		hQ2_vs_x_B_rb.setTitleY("Q^2 [GeV^2]");
		cankin.getPad(0).setTitleFontSize(32);
		cankin.getPad(0).setAxisTitleFontSize(32);
		cankin.getPad(0).setAxisLabelFontSize(24);
		cankin.getPad(0).setStatBoxFontSize(18);
		cankin.draw(hQ2_vs_x_B_rb, "same");
		cankin.cd(1);
		cankin.setFont("Arial");
		hnegt_vs_phi_trento.setTitle("-t vs. #phi");
		hnegt_vs_phi_trento.setTitleX("#phi [#degree]");
		hnegt_vs_phi_trento.setTitleY("-t [GeV^2]");
		cankin.getPad(1).setTitleFontSize(32);
		cankin.getPad(1).setAxisTitleFontSize(32);
		cankin.getPad(1).setAxisLabelFontSize(24);
		cankin.getPad(1).setStatBoxFontSize(18);
		cankin.draw(hnegt_vs_phi_trento, "same");
		
		JFrame framex_B = new JFrame("x_B");
		framex_B.setSize(500, 500);
		EmbeddedCanvas canx_B = new EmbeddedCanvas();
		framex_B.add(canx_B);
		framex_B.setLocationRelativeTo(null);
		framex_B.setVisible(true);
		canx_B.divide(1, 1);
		canx_B.cd(0);
		canx_B.setFont("Arial");
		hx_B_rb.setTitle("x_B");
		hx_B_rb.setTitleX("x_B");
		hx_B_rb.setTitleY("Count");
		hx_B_rb.setOptStat(10);
		canx_B.getPad(0).setTitleFontSize(32);
		canx_B.getPad(0).setAxisTitleFontSize(32);
		canx_B.getPad(0).setAxisLabelFontSize(24);
		canx_B.getPad(0).setStatBoxFontSize(18);
		canx_B.draw(hx_B_rb, "same");
		
		for(int bini = 0; bini < 1000; bini++){
			if(histGroups_x_B_bin.getItem(0).getEntries() < (hx_B_rb.getEntries()*1/8)){
				for(int filli = 0; filli < hx_B_rb.getBinContent(bini); filli++){
					histGroups_x_B_bin.getItem(0).fill(0.0008*bini);
				}
			}
			else if(histGroups_x_B_bin.getItem(1).getEntries() < (2*hx_B_rb.getEntries()*1/8)){
				for(int filli = 0; filli < hx_B_rb.getBinContent(bini); filli++){
					histGroups_x_B_bin.getItem(1).fill(0.0008*bini);
				}
			}
			else if(histGroups_x_B_bin.getItem(2).getEntries() < (3*hx_B_rb.getEntries()*1/8)){
				for(int filli = 0; filli < hx_B_rb.getBinContent(bini); filli++){
					histGroups_x_B_bin.getItem(2).fill(0.0008*bini);
				}
			}
			else if(histGroups_x_B_bin.getItem(3).getEntries() < (2*hx_B_rb.getEntries()*1/8)){
				for(int filli = 0; filli < hx_B_rb.getBinContent(bini); filli++){
					histGroups_x_B_bin.getItem(3).fill(0.0008*bini);
				}
			}
		}
		/*
		for(int bini = 0; bini < 1000; bini++){
			if(histGroups_x_B_bin.getItem(0).getEntries() < (hx_B_rb.getEntries()*4/18)){
				for(int filli = 0; filli < hx_B_rb.getBinContent(bini); filli++){
					histGroups_x_B_bin.getItem(0).fill(0.0008*bini);
				}
			}
			else if(histGroups_x_B_bin.getItem(1).getEntries() < (2*hx_B_rb.getEntries()*3/18)){
				for(int filli = 0; filli < hx_B_rb.getBinContent(bini); filli++){
					histGroups_x_B_bin.getItem(1).fill(0.0008*bini);
				}
			}
			else if(histGroups_x_B_bin.getItem(2).getEntries() < (3*hx_B_rb.getEntries()*2/18)){
				for(int filli = 0; filli < hx_B_rb.getBinContent(bini); filli++){
					histGroups_x_B_bin.getItem(2).fill(0.0008*bini);
				}
			}
			else if(histGroups_x_B_bin.getItem(3).getEntries() < (2*hx_B_rb.getEntries()*1/18)){
				for(int filli = 0; filli < hx_B_rb.getBinContent(bini); filli++){
					histGroups_x_B_bin.getItem(3).fill(0.0008*bini);
				}
			}
		}
		*/
		JFrame frame_x_B_bin = new JFrame("x_B Bins");
		frame_x_B_bin.setSize(1000, 1000);
		EmbeddedCanvas can_x_B_bin = new EmbeddedCanvas();
		frame_x_B_bin.add(can_x_B_bin);
		frame_x_B_bin.setLocationRelativeTo(null);
		frame_x_B_bin.setVisible(true);
		can_x_B_bin.divide(2, 2);
		for(int ih = 0; ih < 4; ih++)
		{
			can_x_B_bin.cd(ih);
			can_x_B_bin.setFont("Arial");
			histGroups_x_B_bin.getItem(ih).setTitle("x_B Bin " + (ih+1));
			histGroups_x_B_bin.getItem(ih).setTitleX("x_B");
			histGroups_x_B_bin.getItem(ih).setOptStat(10);
			if(ih == 0) can_x_B_bin.getPad(ih).getAxisX().setRange(0.1, 0.2);
			if(ih == 1) can_x_B_bin.getPad(ih).getAxisX().setRange(0.15, 0.25);
			if(ih == 2) can_x_B_bin.getPad(ih).getAxisX().setRange(0.2, 0.3);
			can_x_B_bin.getPad(ih).setAxisTitleFontSize(32);
			can_x_B_bin.getPad(ih).setTitleFontSize(32);
			can_x_B_bin.getPad(ih).setAxisLabelFontSize(24);
			can_x_B_bin.getPad(ih).setStatBoxFontSize(18);
			can_x_B_bin.draw(histGroups_x_B_bin.getItem(ih), "same");
		}
		
		JFrame frame_x_B_bin_Q2 = new JFrame("x_B Bins Q^2");
		frame_x_B_bin_Q2.setSize(1000, 1000);
		EmbeddedCanvas can_x_B_bin_Q2 = new EmbeddedCanvas();
		frame_x_B_bin_Q2.add(can_x_B_bin_Q2);
		frame_x_B_bin_Q2.setLocationRelativeTo(null);
		frame_x_B_bin_Q2.setVisible(true);
		can_x_B_bin_Q2.divide(2, 2);
		for(int ih = 0; ih < 4; ih++)
		{
			can_x_B_bin_Q2.cd(ih);
			can_x_B_bin_Q2.setFont("Arial");
			histGroups_x_B_bin_Q2.getItem(ih).setTitle("x_B Bin Q^2 " + (ih+1));
			histGroups_x_B_bin_Q2.getItem(ih).setTitleX("Q^2");
			histGroups_x_B_bin_Q2.getItem(ih).setOptStat(10);
			can_x_B_bin_Q2.getPad(ih).setAxisTitleFontSize(32);
			can_x_B_bin_Q2.getPad(ih).setTitleFontSize(32);
			can_x_B_bin_Q2.getPad(ih).setAxisLabelFontSize(24);
			can_x_B_bin_Q2.getPad(ih).setStatBoxFontSize(18);
			can_x_B_bin_Q2.draw(histGroups_x_B_bin_Q2.getItem(ih), "same");
		}
		
		for(int bini = 0; bini < 1000; bini++){
			if(histGroups_Q2_x_B_bin.getItem(0).getEntries() < histGroups_x_B_bin_Q2.getItem(0).getEntries()/1){
				for(int filli = 0; filli < histGroups_x_B_bin_Q2.getItem(0).getBinContent(bini); filli++){
					histGroups_Q2_x_B_bin.getItem(0).fill(0.008*bini);
				}
			}
		}
		for(int bini = 0; bini < 1000; bini++){
			if(histGroups_Q2_x_B_bin.getItem(1).getEntries() < histGroups_x_B_bin_Q2.getItem(1).getEntries()/2){
				for(int filli = 0; filli < histGroups_x_B_bin_Q2.getItem(1).getBinContent(bini); filli++){
					histGroups_Q2_x_B_bin.getItem(1).fill(0.008*bini);
				}
			}
			else if(histGroups_Q2_x_B_bin.getItem(2).getEntries() < histGroups_x_B_bin_Q2.getItem(1).getEntries()/2){
				for(int filli = 0; filli < histGroups_x_B_bin_Q2.getItem(1).getBinContent(bini); filli++){
					histGroups_Q2_x_B_bin.getItem(2).fill(0.008*bini);
				}
			}
		}
		for(int bini = 0; bini < 1000; bini++){
			if(histGroups_Q2_x_B_bin.getItem(3).getEntries() < histGroups_x_B_bin_Q2.getItem(2).getEntries()/3){
				for(int filli = 0; filli < histGroups_x_B_bin_Q2.getItem(2).getBinContent(bini); filli++){
					histGroups_Q2_x_B_bin.getItem(3).fill(0.008*bini);
				}
			}
			else if(histGroups_Q2_x_B_bin.getItem(4).getEntries() < histGroups_x_B_bin_Q2.getItem(2).getEntries()/3){
				for(int filli = 0; filli < histGroups_x_B_bin_Q2.getItem(2).getBinContent(bini); filli++){
					histGroups_Q2_x_B_bin.getItem(4).fill(0.008*bini);
				}
			}
			else if(histGroups_Q2_x_B_bin.getItem(5).getEntries() < histGroups_x_B_bin_Q2.getItem(2).getEntries()/3){
				for(int filli = 0; filli < histGroups_x_B_bin_Q2.getItem(2).getBinContent(bini); filli++){
					histGroups_Q2_x_B_bin.getItem(5).fill(0.008*bini);
				}
			}
		}
		for(int bini = 0; bini < 1000; bini++){
			if(histGroups_Q2_x_B_bin.getItem(6).getEntries() < histGroups_x_B_bin_Q2.getItem(3).getEntries()/2){
				for(int filli = 0; filli < histGroups_x_B_bin_Q2.getItem(3).getBinContent(bini); filli++){
					histGroups_Q2_x_B_bin.getItem(6).fill(0.008*bini);
				}
			}
			else if(histGroups_Q2_x_B_bin.getItem(7).getEntries() < histGroups_x_B_bin_Q2.getItem(3).getEntries()/2){
				for(int filli = 0; filli < histGroups_x_B_bin_Q2.getItem(3).getBinContent(bini); filli++){
					histGroups_Q2_x_B_bin.getItem(7).fill(0.008*bini);
				}
			}
		}
						
		JFrame frame_Q2_x_B_bin = new JFrame("Q^2-x_B Bin");
		frame_Q2_x_B_bin.setSize(2000, 1000);
		EmbeddedCanvas can_Q2_x_B_bin = new EmbeddedCanvas();
		frame_Q2_x_B_bin.add(can_Q2_x_B_bin);
		frame_Q2_x_B_bin.setLocationRelativeTo(null);
		frame_Q2_x_B_bin.setVisible(true);
		can_Q2_x_B_bin.divide(4, 2);
		for(int ih = 0; ih < 8; ih++)
		{
			can_Q2_x_B_bin.cd(ih);
			can_Q2_x_B_bin.setFont("Arial");
			histGroups_Q2_x_B_bin.getItem(ih).setTitle("Q^2-x_B Bin" + (ih+1));
			histGroups_Q2_x_B_bin.getItem(ih).setTitleX("Q^2");
			histGroups_Q2_x_B_bin.getItem(ih).setOptStat(10);
			can_Q2_x_B_bin.getPad(ih).setAxisTitleFontSize(32);
			can_Q2_x_B_bin.getPad(ih).setTitleFontSize(32);
			can_Q2_x_B_bin.getPad(ih).setAxisLabelFontSize(24);
			can_Q2_x_B_bin.getPad(ih).setStatBoxFontSize(18);
			can_Q2_x_B_bin.draw(histGroups_Q2_x_B_bin.getItem(ih), "same");
		}
		
		JFrame frame_negt_vs_phi_bin = new JFrame("-t Bin");
		frame_negt_vs_phi_bin.setSize(2000, 1000);
		EmbeddedCanvas can_negt_vs_phi_bin = new EmbeddedCanvas();
		frame_negt_vs_phi_bin.add(can_negt_vs_phi_bin);
		frame_negt_vs_phi_bin.setLocationRelativeTo(null);
		frame_negt_vs_phi_bin.setVisible(true);
		can_negt_vs_phi_bin.divide(4, 2);
		for(int ih = 0; ih < 8; ih++)
		{
			can_negt_vs_phi_bin.cd(ih);
			can_negt_vs_phi_bin.setFont("Arial");
			histGroups_negt_vs_phi_trento_bin.getItem(ih).setTitle("-t vs. #phi Bin " + (ih+1));
			histGroups_negt_vs_phi_trento_bin.getItem(ih).setTitleX("#phi [#degree]");
			histGroups_negt_vs_phi_trento_bin.getItem(ih).setTitleY("-t [GeV^2]");
			can_negt_vs_phi_bin.getPad(ih).setAxisTitleFontSize(32);
			can_negt_vs_phi_bin.getPad(ih).setTitleFontSize(32);
			can_negt_vs_phi_bin.getPad(ih).setAxisLabelFontSize(24);
			can_negt_vs_phi_bin.getPad(ih).setStatBoxFontSize(18);
			can_negt_vs_phi_bin.draw(histGroups_negt_vs_phi_trento_bin.getItem(ih), "same");
		}
		
		JFrame framenegt = new JFrame("-t");
		framenegt.setSize(500, 500);
		EmbeddedCanvas cannegt = new EmbeddedCanvas();
		framenegt.add(cannegt);
		framenegt.setLocationRelativeTo(null);
		framenegt.setVisible(true);
		cannegt.divide(1, 1);
		cannegt.cd(0);
		cannegt.setFont("Arial");
		hnegt.setTitle("-t");
		hnegt.setTitleX("-t [GeV^2]");
		hnegt.setTitleY("Count");
		hnegt.setOptStat(10);
		cannegt.getPad(0).setTitleFontSize(32);
		cannegt.getPad(0).setAxisTitleFontSize(32);
		cannegt.getPad(0).setAxisLabelFontSize(24);
		cannegt.getPad(0).setStatBoxFontSize(18);
		cannegt.draw(hnegt, "same");
		
		for(int bini = 0; bini < 1000; bini++){
			if(histGroups_negt_bin.getItem(0).getEntries() < (hnegt.getEntries()/5)){
				for(int filli = 0; filli < hnegt.getBinContent(bini); filli++){
					histGroups_negt_bin.getItem(0).fill(0.003*bini);
				}
			}
			else if(histGroups_negt_bin.getItem(1).getEntries() < (hnegt.getEntries()/5)){
				for(int filli = 0; filli < hnegt.getBinContent(bini); filli++){
					histGroups_negt_bin.getItem(1).fill(0.003*bini);
				}
			}
			else if(histGroups_negt_bin.getItem(2).getEntries() < (hnegt.getEntries()/5)){
				for(int filli = 0; filli < hnegt.getBinContent(bini); filli++){
					histGroups_negt_bin.getItem(2).fill(0.003*bini);
				}
			}
			else if(histGroups_negt_bin.getItem(3).getEntries() < (hnegt.getEntries()/5)){
				for(int filli = 0; filli < hnegt.getBinContent(bini); filli++){
					histGroups_negt_bin.getItem(3).fill(0.003*bini);
				}
			}
			else if(histGroups_negt_bin.getItem(4).getEntries() < (hnegt.getEntries()/5)){
				for(int filli = 0; filli < hnegt.getBinContent(bini); filli++){
					histGroups_negt_bin.getItem(4).fill(0.003*bini);
				}
			}
		}

		JFrame frame_negt_bin = new JFrame("-t Bin");
		frame_negt_bin.setSize(1500, 1000);
		EmbeddedCanvas can_negt_bin = new EmbeddedCanvas();
		frame_negt_bin.add(can_negt_bin);
		frame_negt_bin.setLocationRelativeTo(null);
		frame_negt_bin.setVisible(true);
		can_negt_bin.divide(3, 2);
		for(int ih = 0; ih < 5; ih++)
		{
			can_negt_bin.cd(ih);
			can_negt_bin.setFont("Arial");
			histGroups_negt_bin.getItem(ih).setTitle("-t Bin" + (ih+1));
			histGroups_negt_bin.getItem(ih).setTitleX("-t [GeV^2]");
			histGroups_negt_bin.getItem(ih).setOptStat(110);
			can_negt_bin.getPad(ih).setAxisTitleFontSize(32);
			can_negt_bin.getPad(ih).setTitleFontSize(32);
			can_negt_bin.getPad(ih).setAxisLabelFontSize(24);
			can_negt_bin.getPad(ih).setStatBoxFontSize(18);
			can_negt_bin.draw(histGroups_negt_bin.getItem(ih), "same");
		}		

		JFrame frame_Q2_x_B_bin_negt = new JFrame("Q^2-x_B Bin -t");
		frame_Q2_x_B_bin_negt.setSize(2000, 1000);
		EmbeddedCanvas can_Q2_x_B_bin_negt = new EmbeddedCanvas();
		frame_Q2_x_B_bin_negt.add(can_Q2_x_B_bin_negt);
		frame_Q2_x_B_bin_negt.setLocationRelativeTo(null);
		frame_Q2_x_B_bin_negt.setVisible(true);
		can_Q2_x_B_bin_negt.divide(4, 2);
		for(int ih = 0; ih < 8; ih++)
		{
			can_Q2_x_B_bin_negt.cd(ih);
			can_Q2_x_B_bin_negt.setFont("Arial");
			histGroups_Q2_x_B_bin_negt.getItem(ih).setTitle("Q^2-x_B Bin " + (ih+1) + " -t");
			histGroups_Q2_x_B_bin_negt.getItem(ih).setTitleX("-t [GeV^2]");
			histGroups_Q2_x_B_bin_negt.getItem(ih).setOptStat(10);
			can_Q2_x_B_bin_negt.getPad(ih).setAxisTitleFontSize(32);
			can_Q2_x_B_bin_negt.getPad(ih).setTitleFontSize(32);
			can_Q2_x_B_bin_negt.getPad(ih).setAxisLabelFontSize(24);
			can_Q2_x_B_bin_negt.getPad(ih).setStatBoxFontSize(18);
			can_Q2_x_B_bin_negt.draw(histGroups_Q2_x_B_bin_negt.getItem(ih), "same");
			
			for(int bini = 0; bini < 1000; bini++){
				if(histGroups_Q2_x_B_negt_bin.getItem(ih, 0).getEntries() < (histGroups_Q2_x_B_bin_negt.getItem(ih).getEntries()/5)){
					for(int filli = 0; filli < histGroups_Q2_x_B_bin_negt.getItem(ih).getBinContent(bini); filli++){
						histGroups_Q2_x_B_negt_bin.getItem(ih, 0).fill(0.003*bini);
					}
				}
				else if(histGroups_Q2_x_B_negt_bin.getItem(ih, 1).getEntries() < (histGroups_Q2_x_B_bin_negt.getItem(ih).getEntries()/5)){
					for(int filli = 0; filli < histGroups_Q2_x_B_bin_negt.getItem(ih).getBinContent(bini); filli++){
						histGroups_Q2_x_B_negt_bin.getItem(ih, 1).fill(0.003*bini);
					}
				}
				else if(histGroups_Q2_x_B_negt_bin.getItem(ih, 2).getEntries() < (histGroups_Q2_x_B_bin_negt.getItem(ih).getEntries()/5)){
					for(int filli = 0; filli < histGroups_Q2_x_B_bin_negt.getItem(ih).getBinContent(bini); filli++){
						histGroups_Q2_x_B_negt_bin.getItem(ih, 2).fill(0.003*bini);
					}
				}
				else if(histGroups_Q2_x_B_negt_bin.getItem(ih, 3).getEntries() < (histGroups_Q2_x_B_bin_negt.getItem(ih).getEntries()/5)){
					for(int filli = 0; filli < histGroups_Q2_x_B_bin_negt.getItem(ih).getBinContent(bini); filli++){
						histGroups_Q2_x_B_negt_bin.getItem(ih, 3).fill(0.003*bini);
					}
				}
				else if(histGroups_Q2_x_B_negt_bin.getItem(ih, 4).getEntries() < (histGroups_Q2_x_B_bin_negt.getItem(ih).getEntries()/5)){
					for(int filli = 0; filli < histGroups_Q2_x_B_bin_negt.getItem(ih).getBinContent(bini); filli++){
						histGroups_Q2_x_B_negt_bin.getItem(ih, 4).fill(0.003*bini);
					}
				}
			}
		}
		
		for(int it = 0; it < 8; it++)
		{
			JFrame frame_Q2_x_B_negt_bin = new JFrame("Q^2-x_B Bin " + (it+1) + " -t Bin");
			frame_Q2_x_B_negt_bin.setSize(1500, 1000);
			EmbeddedCanvas can_Q2_x_B_negt_bin = new EmbeddedCanvas();
			frame_Q2_x_B_negt_bin.add(can_Q2_x_B_negt_bin);
			frame_Q2_x_B_negt_bin.setLocationRelativeTo(null);
			frame_Q2_x_B_negt_bin.setVisible(true);
			can_Q2_x_B_negt_bin.divide(3, 2);
			for(int ih = 0; ih < 5; ih++)
			{
				can_Q2_x_B_negt_bin.cd(ih);
				can_Q2_x_B_negt_bin.setFont("Arial");
				histGroups_Q2_x_B_negt_bin.getItem(it, ih).setTitle("Q^2-x_B Bin " + (it+1) + " -t Bin" + (ih+1));
				histGroups_Q2_x_B_negt_bin.getItem(it, ih).setTitleX("-t [GeV^2]");
				histGroups_Q2_x_B_negt_bin.getItem(it, ih).setOptStat(110);
				can_Q2_x_B_negt_bin.getPad(ih).setAxisTitleFontSize(32);
				can_Q2_x_B_negt_bin.getPad(ih).setTitleFontSize(32);
				can_Q2_x_B_negt_bin.getPad(ih).setAxisLabelFontSize(24);
				can_Q2_x_B_negt_bin.getPad(ih).setStatBoxFontSize(18);
				can_Q2_x_B_negt_bin.draw(histGroups_Q2_x_B_negt_bin.getItem(it, ih), "same");
			}
		}
	}
}