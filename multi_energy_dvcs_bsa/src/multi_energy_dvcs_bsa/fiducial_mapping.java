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

public class fiducial_mapping {
	
	static H2F hy_vs_x = new H2F("y_vs_x", "y_vs_x", 500, -500, 500, 500, -500, 500);
	static H2F htheta_vs_phi = new H2F("theta_vs_phi", "theta_vs_phi", 500, -180, 180, 500, 0, 45);
	
	static IndexedList<H2F> histGroups_ecal_det_y_vs_x = new IndexedList<H2F>(1);
	public static void ecal_det_y_vs_x_histos() {
		for(int idet = 0; idet < 3; idet++) 
		{
			H2F ecal_det_y_vs_x = new H2F("ecal_det_y_vs_x", "ecal_det_y_vs_x", 500, -500, 500, 500, -500, 500);
			histGroups_ecal_det_y_vs_x.add(ecal_det_y_vs_x, idet);
		}
	}
	
	static IndexedList<H1F> histGroups_ecal_channel = new IndexedList<H1F>(3);
	public static void ecal_channel_histos() {
		for(int isec = 0; isec  < 6; isec++) {
			for(int idet = 0; idet < 3; idet++) {
				for(int idir = 0; idir < 3; idir++)
				{
					H1F ecal_channel = new H1F("ecal_channel", 900, 0, 450);
					histGroups_ecal_channel.add(ecal_channel, isec, idet, idir);
				}
			}
		}
	}
	
	static H2F hftcal_det_y_vs_x = new H2F("ftcal_det_y_vs_x", "ftcal_det_y_vs_x", 500, -17.5, 17.5, 500, -17.5, 17.5);
	static H2F hftcal_det_theta_vs_phi = new H2F("ftcal_det_theta_vs_phi", "ftcal_det_theta_vs_phi", 500, -180, 180, 500, 0, 7);
	static H2F hftcal_det_cryst_y_vs_x = new H2F("ftcal_det_cryst_y_vs_x", "ftcal_det_cryst_y_vs_x", 24, -18.36, 18.36, 24, -18.36, 18.36);
	
	static HipoDataSource reader = new HipoDataSource();
	
	static double rc = Math.PI/180;
	
	static void processEvent(DataEvent event, int eventCounter)
	{
		if(event.hasBank("REC::Particle"))
		{
			HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
			int gammacount = 0;
			int gamma_pindex = -1;
			double cX = 0;
			double cY = 0;
			boolean fid_ecal_gamma = false;
			boolean fid_ftcal_gamma = false;
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
			LorentzVector p_gamma = new LorentzVector(0, 0, 0, 0);
			for(int i = 0; i < rec.rows(); i++)
			{
				int pid = rec.getInt("pid", i);
				float px = rec.getFloat("px", i);
				float py = rec.getFloat("py", i);
				float pz = rec.getFloat("pz", i);
				float p = (float) Math.sqrt((px*px)+(py*py)+(pz*pz));
				if(pid == 22)
				{				
					p_gamma = new LorentzVector(px, py, pz, p);
					gamma_pindex = i;
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
								if(cX > cX_split[sector-1] && cY < s_left[sector-1]*(cX-t_left[sector-1]) && cY > s_right[sector-1]*(cX-t_right[sector-1]))
									fid_ecal_gamma = true;
								else if(cX < cX_split[sector-1] && cY < q_left[sector-1]*(cX-r_left[sector-1]) && cY > q_right[sector-1]*(cX-r_right[sector-1]))
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
							
							if(fid_ecal_gamma == true)
							{
								if(pindex == gamma_pindex && detector == 7)
								{
									histGroups_ecal_det_y_vs_x.getItem((layer-1)/3).fill(x,y);
									histGroups_ecal_channel.getItem(sector-1, (layer-1)/3, 0).fill(lu);
									histGroups_ecal_channel.getItem(sector-1, (layer-1)/3, 1).fill(lv);
									histGroups_ecal_channel.getItem(sector-1, (layer-1)/3, 2).fill(lw);
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
							if(pindex == gamma_pindex && detector == 10 && layer == 1)
							{
								/*if(ftcal_det.rho() > 8.8 && ftcal_det.rho() > -15.5) fid_ftcal_gamma = true;
								if((x+8.5)*(x+8.5)+(y-10)*(y-10) < 1.5*1.5 
									|| (x+6)*(x+6)/(2.3*2.3)+(y+13)*(y+13)/(1.6*1.6) < 1
									|| (x-4)*(x-4)+(y+6.8)*(y+6.8) < 1.5*1.5) fid_ftcal_gamma = false;*/
								if(ftcal_det.rho() > 8.8 && !(ftcal_det.rho() > 15.5)
									&& (x+8.5)*(x+8.5)+(y-10)*(y-10) > 1.5*1.5 
									&& (x+6)*(x+6)/(2.3*2.3)+(y+13)*(y+13)/(1.6*1.6) > 1
									&& (x-3.8)*(x-3.8)+(y+6.8)*(y+6.8) > 1.7*1.7) fid_ftcal_gamma = true;
							}
						//	if(fid_ftcal_gamma == true)
							{
								if(pindex == gamma_pindex && detector == 10 && layer == 1)
								{
									double ftcal_det_phi = 180*ftcal_det.phi()/Math.PI;
									double ftcal_det_theta = 180*ftcal_det.theta()/Math.PI;
									hftcal_det_y_vs_x.fill(x, y);
									hftcal_det_theta_vs_phi.fill(ftcal_det_phi, ftcal_det_theta);
									hftcal_det_cryst_y_vs_x.fill(x, y);
								}
							}
						}
					}
					if(fid_ecal_gamma == true || fid_ftcal_gamma == true)
					{
						double x_proj = 612.5*p_gamma.px()/p_gamma.pz();
						double y_proj = 612.5*p_gamma.py()/p_gamma.pz();
						double phi = 180*p_gamma.phi()/Math.PI;
						double theta = 180*p_gamma.theta()/Math.PI;
						hy_vs_x.fill(x_proj, y_proj);
						htheta_vs_phi.fill(phi, theta);
					}
				}
			}	
		}
	}
	
	public static void main(String[] args) {
		
		ecal_channel_histos();
		ecal_det_y_vs_x_histos();
		
	//	reader.open("C:/Users/joshtanj/Documents/download/skim_epg_bank_merged_6535MeV_skim16.hipo");
		reader.open("C:/Users/joshtanj/Documents/download/skim_epg_bank_merged_7546MeV_skim16.hipo");
		
		int eventCounter = 1;
		while(reader.hasEvent())// && eventCounter < 20000000)
		{
			processEvent(reader.getNextEvent(), eventCounter);
			if(eventCounter%50000 == 0) System.out.println("Event: " + eventCounter);
			eventCounter++;
		}
		
		JFrame frame = new JFrame("Planar Hits");
		frame.setSize(1500, 750);
		EmbeddedCanvas can = new EmbeddedCanvas();
		frame.add(can);
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		can.divide(2, 1);
		can.cd(0);
		can.setFont("Arial");
		hy_vs_x.setTitle("y vs. x");
		hy_vs_x.setTitleX("x [cm]");
		hy_vs_x.setTitleY("y [cm]");
		can.getPad(0).setTitleFontSize(32);
		can.getPad(0).setAxisTitleFontSize(32);
		can.getPad(0).setAxisLabelFontSize(24);
		can.draw(hy_vs_x, "same");
		can.cd(1);
		can.setFont("Arial");
		htheta_vs_phi.setTitle("#theta vs. #phi");
		htheta_vs_phi.setTitleX("#phi [deg]");
		htheta_vs_phi.setTitleY("#theta [deg]");
		can.getPad(1).setTitleFontSize(32);
		can.getPad(1).setAxisTitleFontSize(32);
		can.getPad(1).setAxisLabelFontSize(24);
		can.draw(htheta_vs_phi, "same");
		
		JFrame frame_ecal_hit = new JFrame("ECAL Hits");
		frame_ecal_hit.setSize(1500, 500);
		EmbeddedCanvas can_ecal_hit = new EmbeddedCanvas();
		frame_ecal_hit.add(can_ecal_hit);
		frame_ecal_hit.setLocationRelativeTo(null);
		frame_ecal_hit.setVisible(true);
		can_ecal_hit.divide(3, 1);
		for(int idet = 0; idet < 3; idet++) 
		{
			can_ecal_hit.cd(idet);
			can_ecal_hit.setFont("Arial");
			histGroups_ecal_det_y_vs_x.getItem(idet).setTitle("ECAL Layer " + (idet*3+1) + " y vs. x");
			histGroups_ecal_det_y_vs_x.getItem(idet).setTitleX("x [cm]");
			histGroups_ecal_det_y_vs_x.getItem(idet).setTitleY("y [cm]");
			can_ecal_hit.getPad(idet).setTitleFontSize(32);
			can_ecal_hit.getPad(idet).setAxisTitleFontSize(32);
			can_ecal_hit.getPad(idet).setAxisLabelFontSize(24);
			can_ecal_hit.draw(histGroups_ecal_det_y_vs_x.getItem(idet), "same");
		}
		
		for(int idet = 0; idet < 3; idet++) {
			for(int idir = 0; idir < 3; idir++) {
				JFrame frame_ecal_channel = new JFrame("ECAL Channels");
				frame_ecal_channel.setSize(1500, 1000);
				EmbeddedCanvas can_ecal_channel = new EmbeddedCanvas();
				frame_ecal_channel.add(can_ecal_channel);
				frame_ecal_channel.setLocationRelativeTo(null);
				frame_ecal_channel.setVisible(true);
				can_ecal_channel.divide(1, 6);
				can_ecal_channel.cd(0);
				can_ecal_channel.setFont("Arial");
				histGroups_ecal_channel.getItem(0, idet, idir).setTitleY("[Sec. 1]");
				histGroups_ecal_channel.getItem(0, idet, idir).setFillColor(32);
				can_ecal_channel.getPad(0).setAxisTitleFontSize(24);
				can_ecal_channel.getPad(0).setAxisLabelFontSize(16);
				can_ecal_channel.getPad(0).getAxisY().setLog(true);
				can_ecal_channel.draw(histGroups_ecal_channel.getItem(0, idet, idir), "same");
				can_ecal_channel.cd(1);
				can_ecal_channel.setFont("Arial");
				histGroups_ecal_channel.getItem(1, idet, idir).setTitleY("[Sec. 2]");
				histGroups_ecal_channel.getItem(1, idet, idir).setFillColor(32);
				can_ecal_channel.getPad(1).setAxisTitleFontSize(24);
				can_ecal_channel.getPad(1).setAxisLabelFontSize(16);
				can_ecal_channel.getPad(1).getAxisY().setLog(true);
				can_ecal_channel.draw(histGroups_ecal_channel.getItem(1, idet, idir), "same");
				can_ecal_channel.cd(2);
				can_ecal_channel.setFont("Arial");
				histGroups_ecal_channel.getItem(2, idet, idir).setTitleY("[Sec. 3]");
				histGroups_ecal_channel.getItem(2, idet, idir).setFillColor(32);
				can_ecal_channel.getPad(2).setAxisTitleFontSize(24);
				can_ecal_channel.getPad(2).setAxisLabelFontSize(16);
				can_ecal_channel.getPad(2).getAxisY().setLog(true);
				can_ecal_channel.draw(histGroups_ecal_channel.getItem(2, idet, idir), "same");
				can_ecal_channel.cd(3);
				can_ecal_channel.setFont("Arial");
				histGroups_ecal_channel.getItem(3, idet, idir).setTitleY("[Sec. 4]");
				histGroups_ecal_channel.getItem(3, idet, idir).setFillColor(32);
				can_ecal_channel.getPad(3).setAxisTitleFontSize(24);
				can_ecal_channel.getPad(3).setAxisLabelFontSize(16);
				can_ecal_channel.getPad(3).getAxisY().setLog(true);
				can_ecal_channel.draw(histGroups_ecal_channel.getItem(3, idet, idir), "same");
				can_ecal_channel.cd(4);
				can_ecal_channel.setFont("Arial");
				histGroups_ecal_channel.getItem(4, idet, idir).setTitleY("[Sec. 5]");
				histGroups_ecal_channel.getItem(4, idet, idir).setFillColor(32);
				can_ecal_channel.getPad(4).setAxisTitleFontSize(24);
				can_ecal_channel.getPad(4).setAxisLabelFontSize(16);
				can_ecal_channel.getPad(4).getAxisY().setLog(true);
				can_ecal_channel.draw(histGroups_ecal_channel.getItem(4, idet, idir), "same");
				can_ecal_channel.cd(5);
				can_ecal_channel.setFont("Arial");
				if(idir == 0) histGroups_ecal_channel.getItem(5, idet, idir).setTitleX("ECAL Layer " + (idet*3+1) + " lu [cm]");
				if(idir == 1) histGroups_ecal_channel.getItem(5, idet, idir).setTitleX("ECAL Layer " + (idet*3+1) + " lv [cm]");
				if(idir == 2) histGroups_ecal_channel.getItem(5, idet, idir).setTitleX("ECAL Layer " + (idet*3+1) + " lw [cm]");
				histGroups_ecal_channel.getItem(5, idet, idir).setTitleY("[Sec. 6]");
				histGroups_ecal_channel.getItem(5, idet, idir).setFillColor(32);
				can_ecal_channel.getPad(5).setAxisTitleFontSize(24);
				can_ecal_channel.getPad(5).setAxisLabelFontSize(16);
				can_ecal_channel.getPad(5).getAxisY().setLog(true);
				can_ecal_channel.draw(histGroups_ecal_channel.getItem(5, idet, idir), "same");
			}
		}
		
		JFrame frame_ftcal_hit = new JFrame("FTCAL Hits");
		frame_ftcal_hit.setSize(1500, 750);
		EmbeddedCanvas can_ftcal_hit = new EmbeddedCanvas();
		frame_ftcal_hit.add(can_ftcal_hit);
		frame_ftcal_hit.setLocationRelativeTo(null);
		frame_ftcal_hit.setVisible(true);
		can_ftcal_hit.divide(2, 1);
		can_ftcal_hit.cd(0);
		can_ftcal_hit.setFont("Arial");
		hftcal_det_y_vs_x.setTitle("FTCAL y vs. x");
		hftcal_det_y_vs_x.setTitleX("x [cm]");
		hftcal_det_y_vs_x.setTitleY("y [cm]");
		can_ftcal_hit.getPad(0).setTitleFontSize(32);
		can_ftcal_hit.getPad(0).setAxisTitleFontSize(32);
		can_ftcal_hit.getPad(0).setAxisLabelFontSize(24);
		can_ftcal_hit.draw(hftcal_det_y_vs_x, "same");
		can_ftcal_hit.cd(1);
		can_ftcal_hit.setFont("Arial");
		hftcal_det_theta_vs_phi.setTitle("FTCAL #theta vs. #phi");
		hftcal_det_theta_vs_phi.setTitleX("#phi [deg]");
		hftcal_det_theta_vs_phi.setTitleY("#theta [deg]");
		can_ftcal_hit.getPad(1).setTitleFontSize(32);
		can_ftcal_hit.getPad(1).setAxisTitleFontSize(32);
		can_ftcal_hit.getPad(1).setAxisLabelFontSize(24);
		can_ftcal_hit.draw(hftcal_det_theta_vs_phi, "same");
		
		JFrame frame_ftcal_cryst_hit = new JFrame("FTCAL Crystal Hits");
		frame_ftcal_cryst_hit.setSize(750, 750);
		EmbeddedCanvas can_ftcal_cryst_hit = new EmbeddedCanvas();
		frame_ftcal_cryst_hit.add(can_ftcal_cryst_hit);
		frame_ftcal_cryst_hit.setLocationRelativeTo(null);
		frame_ftcal_cryst_hit.setVisible(true);
		can_ftcal_cryst_hit.divide(1, 1);
		can_ftcal_cryst_hit.cd(0);
		can_ftcal_cryst_hit.setFont("Arial");
		hftcal_det_cryst_y_vs_x.setTitle("FTCAL Crystal y vs. x");
		hftcal_det_cryst_y_vs_x.setTitleX("x [cm]");
		hftcal_det_cryst_y_vs_x.setTitleY("y [cm]");
		can_ftcal_cryst_hit.getPad(0).setTitleFontSize(32);
		can_ftcal_cryst_hit.getPad(0).setAxisTitleFontSize(32);
		can_ftcal_cryst_hit.getPad(0).setAxisLabelFontSize(24);
		can_ftcal_cryst_hit.draw(hftcal_det_cryst_y_vs_x, "same");
	}
}