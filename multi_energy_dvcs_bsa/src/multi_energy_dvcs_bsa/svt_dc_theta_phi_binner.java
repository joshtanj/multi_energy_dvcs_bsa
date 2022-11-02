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

public class svt_dc_theta_phi_binner {
	
	static HipoDataSource reader = new HipoDataSource();
	static H2F hdc_det_y_vs_x = new H2F("dc_det_y_vs_x", "dc_det_y_vs_x", 500, -200, 200, 500, -200, 200);
	static H1F hdc_det_theta = new H1F("dc_det_theta", "dc_det_theta", 200, 17, 47);
	static H2F hdc_det_theta_vs_phi = new H2F("dc_det_theta_vs_phi", "dc_det_theta_vs_phi", 1500, -150, 210, 500, 17, 47);
	static H2F hdc_det_theta_vs_phi_rot = new H2F("dc_det_theta_vs_phi_rot", "dc_det_theta_vs_phi_rot", 500, -60, 60, 500, 17, 47);
	static H2F hdc_det_theta_vs_phi_rot_mir = new H2F("dc_det_theta_vs_phi_rot_mir", "dc_det_theta_vs_phi_rot_mir", 500, -60, 60, 500, 17, 47);
	static H1F hdc_det_phi_rot_mir_prj = new H1F("dc_det_phi_rot_mir_prj", "dc_det_phi_rot_mir_prj", 200, -30, 30);
	
	static H2F hsvt_det_y_vs_x = new H2F("svt_det_y_vs_x", "svt_det_y_vs_x", 500, -50, 50, 500, -50, 50);
	static H1F hsvt_det_theta = new H1F("svt_det_theta", "svt_det_theta", 200, 32, 87);
	static H2F hsvt_det_theta_vs_phi = new H2F("svt_det_theta_vs_phi", "svt_det_theta_vs_phi", 1500, -90, 270, 500, 32, 87);
	static H2F hsvt_det_theta_vs_phi_rot = new H2F("svt_det_theta_vs_phi_rot", "svt_det_theta_vs_phi_rot", 500, -120, 120, 500, 32, 87);
	static H1F hsvt_det_phi_rot_prj = new H1F("svt_det_phi_rot_prj", "svt_det_phi_rot_prj", 200, -120, 120);
	
	static IndexedList<H1F> histGroups_dc_det_theta_bin= new IndexedList<H1F>(1);
	public static void dc_det_theta_bin_histos() {
		for(int ibin = 0; ibin < 4; ibin++) {
			H1F dc_det_theta_bin = new H1F("dc_det_theta_bin", "dc_det_theta_bin", 240, 17, 42);
			histGroups_dc_det_theta_bin.add(dc_det_theta_bin, ibin);
		}
	}
	
	static H1F hdc_det_phi_rot_mir_prj_bin1 = new H1F("dc_det_phi_rot_mir_prj_bin1", "dc_det_phi_rot_mir_prj_bin1", 200, -30, 30);
	static H1F hdc_det_phi_rot_mir_prj_bin2 = new H1F("dc_det_phi_rot_mir_prj_bin2", "dc_det_phi_rot_mir_prj_bin2", 200, -30, 30);
	
	static IndexedList<H1F> histGroups_svt_det_theta_bin= new IndexedList<H1F>(1);
	public static void svt_det_theta_bin_histos() {
		for(int ibin = 0; ibin < 9; ibin++) {
			H1F svt_det_theta_bin = new H1F("svt_det_theta_bin", "svt_det_theta_bin", 240, 32, 87);
			histGroups_svt_det_theta_bin.add(svt_det_theta_bin, ibin);
		}
	}
	
	static H1F hsvt_det_phi_rot_prj_bin1 = new H1F("svt_det_phi_rot_mir_prj_bin1", "svt_det_phi_rot_mir_prj_bin1", 200, -60, 60);
	static H1F hsvt_det_phi_rot_prj_bin2 = new H1F("svt_det_phi_rot_mir_prj_bin2", "svt_det_phi_rot_mir_prj_bin2", 200, -60, 60);
	static H1F hsvt_det_phi_rot_prj_bin3 = new H1F("svt_det_phi_rot_mir_prj_bin3", "svt_det_phi_rot_mir_prj_bin3", 200, -60, 60);
	static H1F hsvt_det_phi_rot_prj_bin4 = new H1F("svt_det_phi_rot_mir_prj_bin4", "svt_det_phi_rot_mir_prj_bin4", 200, -60, 60);
	
	static double dtr = Math.PI/180;
	
	static void processEvent(DataEvent event, int eventCounter)
	{
		if(event.hasBank("REC::Particle") && event.hasBank("REC::Track") && event.hasBank("REC::Traj"))
		{
			HipoDataBank rec = (HipoDataBank) event.getBank("REC::Particle");
			int ecount = 0;
			int protoncount = 0;
			int e_pindex = -1;
			int proton_pindex = -1;
			float theta_deg_proton = 0;
			float phi_deg_proton = -200;
			float phi_rot_deg_proton = -200;
			double[] phi_rot  = new double[] {0, 60, 120, 180, 240, 300};
			double[] phi_crot  = new double[] {90, 210, 330};
			LorentzVector p_e = new LorentzVector(0, 0, 0, 0.000511);
			LorentzVector p_proton = new LorentzVector(0, 0, 0, 0.938272);
		//	double E = 6.535;
			double E = 7.546;
			Vector3D v_e = new Vector3D (0, 0, 0);
			Vector3D v_proton = new Vector3D (0, 0, 0);
			float cdchi2_proton = 100;
			short cdNDF_proton = 100;
			Vector3D dc_hit_proton = new Vector3D (0, 0, 0);
			Vector3D dc_hit_rot_proton = new Vector3D (0, 0, 0);
			Vector3D svt_hit_proton = new Vector3D (0, 0, 0);
			Vector3D svt_hit_rot_proton = new Vector3D (0, 0, 0);
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
			byte proton_detector = 0;
			byte proton_sector = 0;
			for(int j = 0; j < rectrac.rows(); j++)
			{
				short pindex = rectrac.getShort("pindex", j);
				byte detector = rectrac.getByte("detector", j);
				byte sector = rectrac.getByte("sector", j);
				float chi2 = rectrac.getFloat("chi2", j);
				short NDF = rectrac.getShort("NDF", j);
				
				if(pindex == proton_pindex)
				{
					proton_detector = detector;
					proton_sector = sector;
					cdchi2_proton = chi2;
					cdNDF_proton = NDF;
				}
			}
			int proton_csector = 0;
			HipoDataBank rectraj = (HipoDataBank) event.getBank("REC::Traj");
			for(int k = 0; k < rectraj.rows(); k++)
			{
				short pindex = rectraj.getShort("pindex", k);
				byte detector = rectraj.getByte("detector", k);
				byte layer = rectraj.getByte("layer", k);
				float x = rectraj.getFloat("x", k);
				float y = rectraj.getFloat("y", k);
				float z = rectraj.getFloat("z", k);
				if(pindex == proton_pindex && detector == 6 && proton_detector == 6 && layer ==  6)
				{
					dc_hit_proton = new Vector3D (x, y, z);
					theta_deg_proton = (float) (57.3*dc_hit_proton.theta());
					phi_deg_proton = (float) (57.3*dc_hit_proton.phi());
					if(phi_deg_proton <= -150) phi_deg_proton = phi_deg_proton+360;
					dc_hit_rot_proton = new Vector3D (x*Math.cos(phi_rot[proton_sector-1]/57.3)+y*Math.sin(phi_rot[proton_sector-1]/57.3),
													y*Math.cos(phi_rot[proton_sector-1]/57.3)-x*Math.sin(phi_rot[proton_sector-1]/57.3), z);
					phi_rot_deg_proton = (float) (57.3*dc_hit_rot_proton.phi());
				}
				else if(pindex == proton_pindex && detector == 5 && proton_detector == 5 && layer == 12)
				{
					svt_hit_proton = new Vector3D (x, y, z);
					theta_deg_proton = (float) (57.3*svt_hit_proton.theta());
					phi_deg_proton = (float) (57.3*svt_hit_proton.phi());
					if(phi_deg_proton <= -90) phi_deg_proton = phi_deg_proton+360;
					if(-90 <= phi_deg_proton && phi_deg_proton < 30) proton_csector = 3;
					else if(30 <= phi_deg_proton && phi_deg_proton < 150) proton_csector = 1;
					else if(150 <= phi_deg_proton && phi_deg_proton < 270) proton_csector = 2;
					svt_hit_rot_proton = new Vector3D (x*Math.cos(phi_crot[proton_csector-1]/57.3)+y*Math.sin(phi_crot[proton_csector-1]/57.3),
													y*Math.cos(phi_crot[proton_csector-1]/57.3)-x*Math.sin(phi_crot[proton_csector-1]/57.3), z);
					phi_rot_deg_proton = (float) (57.3*svt_hit_rot_proton.phi());
				}
			}
			
			if(ecount == 1 && protoncount == 1
					&& v_e.z() > -12 && v_e.z() < 7 && Math.abs(v_e.z()-v_proton.z()) < 2.5 + (2.5/p_proton.p())
					&& ((proton_detector == 6) || (proton_detector == 5 && cdNDF_proton < 10 && cdchi2_proton/cdNDF_proton < 30)))
				//	&& 57.3*p_proton.theta() < 75)
			{
				if(proton_detector == 6)
				{
					hdc_det_y_vs_x.fill(dc_hit_proton.x(), dc_hit_proton.y());
					hdc_det_theta.fill(theta_deg_proton);
				//	if(theta_deg_proton > 7 && theta_deg_proton < 33)
					{
						hdc_det_theta_vs_phi.fill(phi_deg_proton, theta_deg_proton);
						hdc_det_theta_vs_phi_rot.fill(phi_rot_deg_proton, theta_deg_proton);
						hdc_det_theta_vs_phi_rot_mir.fill(Math.abs(phi_rot_deg_proton), theta_deg_proton);
						hdc_det_phi_rot_mir_prj.fill(Math.abs(phi_rot_deg_proton));
					}
					double[] theta_bnd  = new double[] {17, 32, 35, 38, 42};
					for(int thbin = 0; thbin < 4; thbin++)
					{
						if(theta_deg_proton > theta_bnd[thbin] && theta_deg_proton <= theta_bnd[thbin+1])
						{
							histGroups_dc_det_theta_bin.getItem(thbin).fill(theta_deg_proton);
							break;
						}
					}
				}
				else if(proton_detector == 5)
				{
					hsvt_det_y_vs_x.fill(svt_hit_proton.x(), svt_hit_proton.y());
					hsvt_det_theta.fill(theta_deg_proton);
				//	if(theta_deg_proton > 7 && theta_deg_proton < 33)
					{
						hsvt_det_theta_vs_phi.fill(phi_deg_proton, theta_deg_proton);
						hsvt_det_theta_vs_phi_rot.fill(phi_rot_deg_proton, theta_deg_proton);
						hsvt_det_phi_rot_prj.fill(phi_rot_deg_proton);
					}
					double[] theta_cbnd  = new double[] {40, 59, 63, 66, 68, 70, 72, 74, 77, 80};
					for(int thbin = 0; thbin < 9; thbin++)
					{
						if(theta_deg_proton > theta_cbnd[thbin] && theta_deg_proton <= theta_cbnd[thbin+1])
						{
							histGroups_svt_det_theta_bin.getItem(thbin).fill(theta_deg_proton);
							break;
						}
					}
					
					if(phi_rot_deg_proton <= -15) hsvt_det_phi_rot_prj_bin1.fill(phi_rot_deg_proton);
					else if(phi_rot_deg_proton > -15 && phi_rot_deg_proton <= 10) hsvt_det_phi_rot_prj_bin2.fill(phi_rot_deg_proton);
					else if(phi_rot_deg_proton > 10 && phi_rot_deg_proton <= 43) hsvt_det_phi_rot_prj_bin3.fill(phi_rot_deg_proton);
					else if(phi_rot_deg_proton > 43) hsvt_det_phi_rot_prj_bin4.fill(phi_rot_deg_proton);
				}
			}
		}
	}
	
	public static void main(String[] args) {
		
		dc_det_theta_bin_histos();
		svt_det_theta_bin_histos();
	
	//	reader.open("C:/Users/joshtanj/Documents/download/merged_skim_e_bank_6535MeV_skim_elastic_5950_5899_5963_5958_5885.hipo");
	//	reader.open("C:/Users/joshtanj/Documents/download/merged_skim_ep_bank_7546MeV_skim_elastic.hipo");
		reader.open("C:/Users/joshtanj/Documents/download/skim_6bank_epg_merged_7546MeV_skim16.hipo");
		
		int eventCounter = 1;
		while(reader.hasEvent())// && eventCounter < 10000000)
		{
			processEvent(reader.getNextEvent(), eventCounter);
			if(eventCounter%50000 == 0) System.out.println("Event: " + eventCounter);
			eventCounter++;
		}
		int Nevent = eventCounter-1;
		System.out.println("Number of events: " + Nevent);
		
		JFrame frame_dc_hit = new JFrame("DC Hits");
		frame_dc_hit.setSize(500, 500);
		EmbeddedCanvas can_dc_hit = new EmbeddedCanvas();
		frame_dc_hit.add(can_dc_hit);
		frame_dc_hit.setLocationRelativeTo(null);
		frame_dc_hit.setVisible(true);
		can_dc_hit.divide(1, 1);
		can_dc_hit.cd(0);
		can_dc_hit.setFont("Arial");
		hdc_det_y_vs_x.setTitle("DC y vs. x");
		hdc_det_y_vs_x.setTitleX("x [cm]");
		hdc_det_y_vs_x.setTitleY("y [cm]");
		can_dc_hit.getPad(0).setTitleFontSize(32);
		can_dc_hit.getPad(0).setAxisTitleFontSize(32);
		can_dc_hit.getPad(0).setAxisLabelFontSize(24);
		can_dc_hit.getPad(0).setStatBoxFontSize(18);
		can_dc_hit.draw(hdc_det_y_vs_x, "same");
		
		JFrame frame_dc_hit_theta = new JFrame("DC Hits #theta");
		frame_dc_hit_theta.setSize(500, 500);
		EmbeddedCanvas can_dc_hit_theta = new EmbeddedCanvas();
		frame_dc_hit_theta.add(can_dc_hit_theta);
		frame_dc_hit_theta.setLocationRelativeTo(null);
		frame_dc_hit_theta.setVisible(true);
		can_dc_hit_theta.divide(1, 1);
		can_dc_hit_theta.cd(0);
		can_dc_hit_theta.setFont("Arial");
		hdc_det_theta.setTitle("DC #theta");
		hdc_det_theta.setTitleX("#theta [#degree]");
		hdc_det_theta.setTitleY("Counts");
		hdc_det_theta.setOptStat(10);
		can_dc_hit_theta.getPad(0).setTitleFontSize(32);
		can_dc_hit_theta.getPad(0).setAxisTitleFontSize(32);
		can_dc_hit_theta.getPad(0).setAxisLabelFontSize(24);
		can_dc_hit_theta.getPad(0).setStatBoxFontSize(18);
		can_dc_hit_theta.draw(hdc_det_theta, "same");
		
		JFrame frame_dc_hit_theta_vs_phi = new JFrame("DC Hits #theta vs. #phi");
		frame_dc_hit_theta_vs_phi.setSize(1500, 500);
		EmbeddedCanvas can_dc_hit_theta_vs_phi = new EmbeddedCanvas();
		frame_dc_hit_theta_vs_phi.add(can_dc_hit_theta_vs_phi);
		frame_dc_hit_theta_vs_phi.setLocationRelativeTo(null);
		frame_dc_hit_theta_vs_phi.setVisible(true);
		can_dc_hit_theta_vs_phi.divide(1, 1);
		can_dc_hit_theta_vs_phi.cd(0);
		can_dc_hit_theta_vs_phi.setFont("Arial");
		hdc_det_theta_vs_phi.setTitle("DC #theta vs #phi");
		hdc_det_theta_vs_phi.setTitleX("#phi [#degree]");
		hdc_det_theta_vs_phi.setTitleY("#theta [#degree]");
		can_dc_hit_theta_vs_phi.getPad(0).setTitleFontSize(32);
		can_dc_hit_theta_vs_phi.getPad(0).setAxisTitleFontSize(32);
		can_dc_hit_theta_vs_phi.getPad(0).setAxisLabelFontSize(24);
		can_dc_hit_theta_vs_phi.getPad(0).setStatBoxFontSize(18);
		can_dc_hit_theta_vs_phi.draw(hdc_det_theta_vs_phi, "same");
		
		JFrame frame_dc_hit_theta_vs_phi_prj = new JFrame("DC Hits #phi Projections");
		frame_dc_hit_theta_vs_phi_prj.setSize(1500, 500);
		EmbeddedCanvas can_dc_hit_theta_vs_phi_prj = new EmbeddedCanvas();
		frame_dc_hit_theta_vs_phi_prj.add(can_dc_hit_theta_vs_phi_prj);
		frame_dc_hit_theta_vs_phi_prj.setLocationRelativeTo(null);
		frame_dc_hit_theta_vs_phi_prj.setVisible(true);
		can_dc_hit_theta_vs_phi_prj.divide(3, 1);
		can_dc_hit_theta_vs_phi_prj.cd(0);
		can_dc_hit_theta_vs_phi_prj.setFont("Arial");
		hdc_det_theta_vs_phi_rot.setTitle("Superpositioned DC #theta vs. #phi");
		hdc_det_theta_vs_phi_rot.setTitleX("#phi [#degree]");
		hdc_det_theta_vs_phi_rot.setTitleY("#theta [#degree]");
		can_dc_hit_theta_vs_phi_prj.getPad(0).setTitleFontSize(32);
		can_dc_hit_theta_vs_phi_prj.getPad(0).setAxisTitleFontSize(32);
		can_dc_hit_theta_vs_phi_prj.getPad(0).setAxisLabelFontSize(24);
		can_dc_hit_theta_vs_phi_prj.getPad(0).setStatBoxFontSize(18);
		can_dc_hit_theta_vs_phi_prj.draw(hdc_det_theta_vs_phi_rot, "same");
		can_dc_hit_theta_vs_phi_prj.cd(1);
		can_dc_hit_theta_vs_phi_prj.setFont("Arial");
		hdc_det_theta_vs_phi_rot_mir.setTitle("Superpositioned Mirrored DC #theta vs. #phi");
		hdc_det_theta_vs_phi_rot_mir.setTitleX("#phi [#degree]");
		hdc_det_theta_vs_phi_rot_mir.setTitleY("#theta [#degree]");
		can_dc_hit_theta_vs_phi_prj.getPad(1).setTitleFontSize(32);
		can_dc_hit_theta_vs_phi_prj.getPad(1).setAxisTitleFontSize(32);
		can_dc_hit_theta_vs_phi_prj.getPad(1).setAxisLabelFontSize(24);
		can_dc_hit_theta_vs_phi_prj.getPad(1).setStatBoxFontSize(18);
		can_dc_hit_theta_vs_phi_prj.draw(hdc_det_theta_vs_phi_rot_mir, "same");
		can_dc_hit_theta_vs_phi_prj.cd(2);
		can_dc_hit_theta_vs_phi_prj.setFont("Arial");
		hdc_det_phi_rot_mir_prj.setTitle("Superpositioned Mirrored DC #phi");
		hdc_det_phi_rot_mir_prj.setTitleX("#phi [#degree]");
		hdc_det_phi_rot_mir_prj.setTitleY("Counts");
		hdc_det_phi_rot_mir_prj.setOptStat(10);
		can_dc_hit_theta_vs_phi_prj.getPad(2).setTitleFontSize(32);
		can_dc_hit_theta_vs_phi_prj.getPad(2).setAxisTitleFontSize(32);
		can_dc_hit_theta_vs_phi_prj.getPad(2).setAxisLabelFontSize(24);
		can_dc_hit_theta_vs_phi_prj.getPad(2).setStatBoxFontSize(18);
		can_dc_hit_theta_vs_phi_prj.draw(hdc_det_phi_rot_mir_prj, "same");
		
		double[] theta_bnd  = new double[] {17, 32, 35, 38, 42};
		
		JFrame frame_dc_hit_theta_bin = new JFrame("DC Hits #theta Bins");
		frame_dc_hit_theta_bin.setSize(1000, 1000);
		EmbeddedCanvas can_dc_hit_theta_bin = new EmbeddedCanvas();
		frame_dc_hit_theta_bin.add(can_dc_hit_theta_bin);
		frame_dc_hit_theta_bin.setLocationRelativeTo(null);
		frame_dc_hit_theta_bin.setVisible(true);
		can_dc_hit_theta_bin.divide(2, 2);	
		for(int thbin = 0; thbin < 4; thbin++)
		{
			can_dc_hit_theta_bin.cd(thbin);
			can_dc_hit_theta_bin.setFont("Arial");
			histGroups_dc_det_theta_bin.getItem(thbin).setTitle("DC " + theta_bnd[thbin] + "#degree < #theta <= "
																	+ theta_bnd[thbin+1] + "#degree");
			histGroups_dc_det_theta_bin.getItem(thbin).setTitleX("#theta [#degree]");
			histGroups_dc_det_theta_bin.getItem(thbin).setTitleY("Counts");
			histGroups_dc_det_theta_bin.getItem(thbin).setOptStat(110);
			can_dc_hit_theta_bin.getPad(thbin).setTitleFontSize(32);
			can_dc_hit_theta_bin.getPad(thbin).setAxisTitleFontSize(32);
			can_dc_hit_theta_bin.getPad(thbin).setAxisLabelFontSize(24);
			can_dc_hit_theta_bin.getPad(thbin).setStatBoxFontSize(18);
			can_dc_hit_theta_bin.draw(histGroups_dc_det_theta_bin.getItem(thbin), "same");
		}
		
		for(int bini = 0; bini < 200; bini++){
			if(hdc_det_phi_rot_mir_prj_bin1.getEntries() < hdc_det_phi_rot_mir_prj.getEntries()/2){
				for(int filli = 0; filli < hdc_det_phi_rot_mir_prj.getBinContent(bini); filli++){
					hdc_det_phi_rot_mir_prj_bin1.fill(-30+0.3*bini);
				}
			}
			else{
				for(int filli = 0; filli < hdc_det_phi_rot_mir_prj.getBinContent(bini); filli++){
					hdc_det_phi_rot_mir_prj_bin2.fill(-30+0.3*bini);
				}
			}
		}
		
		JFrame frame_dc_hit_phi_bin = new JFrame("DC Hits #phi Bins");
		frame_dc_hit_phi_bin.setSize(1000, 500);
		EmbeddedCanvas can_dc_hit_phi_bin = new EmbeddedCanvas();
		frame_dc_hit_phi_bin.add(can_dc_hit_phi_bin);
		frame_dc_hit_phi_bin.setLocationRelativeTo(null);
		frame_dc_hit_phi_bin.setVisible(true);
		can_dc_hit_phi_bin.divide(2, 1);
		can_dc_hit_phi_bin.cd(0);
		can_dc_hit_phi_bin.setFont("Arial");
		hdc_det_phi_rot_mir_prj_bin1.setTitle("DC #phi Bin 1");
		hdc_det_phi_rot_mir_prj_bin1.setTitleX("#phi [#degree]");
		hdc_det_phi_rot_mir_prj_bin1.setTitleY("Counts");
		hdc_det_phi_rot_mir_prj_bin1.setOptStat(110);
		can_dc_hit_phi_bin.getPad(0).setTitleFontSize(32);
		can_dc_hit_phi_bin.getPad(0).setAxisTitleFontSize(32);
		can_dc_hit_phi_bin.getPad(0).setAxisLabelFontSize(24);
		can_dc_hit_phi_bin.getPad(0).setStatBoxFontSize(18);
		can_dc_hit_phi_bin.draw(hdc_det_phi_rot_mir_prj_bin1, "same");
		can_dc_hit_phi_bin.cd(1);
		can_dc_hit_phi_bin.setFont("Arial");
		hdc_det_phi_rot_mir_prj_bin2.setTitle("DC #phi Bin 2");
		hdc_det_phi_rot_mir_prj_bin2.setTitleX("#phi [#degree]");
		hdc_det_phi_rot_mir_prj_bin2.setTitleY("Counts");
		hdc_det_phi_rot_mir_prj_bin2.setOptStat(110);
		can_dc_hit_phi_bin.getPad(1).setTitleFontSize(32);
		can_dc_hit_phi_bin.getPad(1).setAxisTitleFontSize(32);
		can_dc_hit_phi_bin.getPad(1).setAxisLabelFontSize(24);
		can_dc_hit_phi_bin.getPad(1).setStatBoxFontSize(18);
		can_dc_hit_phi_bin.draw(hdc_det_phi_rot_mir_prj_bin2, "same");
		
		JFrame frame_svt_hit = new JFrame("SVT Hits");
		frame_svt_hit.setSize(500, 500);
		EmbeddedCanvas can_svt_hit = new EmbeddedCanvas();
		frame_svt_hit.add(can_svt_hit);
		frame_svt_hit.setLocationRelativeTo(null);
		frame_svt_hit.setVisible(true);
		can_svt_hit.divide(1, 1);
		can_svt_hit.cd(0);
		can_svt_hit.setFont("Arial");
		hsvt_det_y_vs_x.setTitle("SVT y vs. x");
		hsvt_det_y_vs_x.setTitleX("x [cm]");
		hsvt_det_y_vs_x.setTitleY("y [cm]");
		can_svt_hit.getPad(0).setTitleFontSize(32);
		can_svt_hit.getPad(0).setAxisTitleFontSize(32);
		can_svt_hit.getPad(0).setAxisLabelFontSize(24);
		can_svt_hit.getPad(0).setStatBoxFontSize(18);
		can_svt_hit.draw(hsvt_det_y_vs_x, "same");
		
		JFrame frame_svt_hit_theta = new JFrame("SVT Hits #theta");
		frame_svt_hit_theta.setSize(500, 500);
		EmbeddedCanvas can_svt_hit_theta = new EmbeddedCanvas();
		frame_svt_hit_theta.add(can_svt_hit_theta);
		frame_svt_hit_theta.setLocationRelativeTo(null);
		frame_svt_hit_theta.setVisible(true);
		can_svt_hit_theta.divide(1, 1);
		can_svt_hit_theta.cd(0);
		can_svt_hit_theta.setFont("Arial");
		hsvt_det_theta.setTitle("SVT #theta");
		hsvt_det_theta.setTitleX("#theta [#degree]");
		hsvt_det_theta.setTitleY("Counts");
		hsvt_det_theta.setOptStat(10);
		can_svt_hit_theta.getPad(0).setTitleFontSize(32);
		can_svt_hit_theta.getPad(0).setAxisTitleFontSize(32);
		can_svt_hit_theta.getPad(0).setAxisLabelFontSize(24);
		can_svt_hit_theta.getPad(0).setStatBoxFontSize(18);
		can_svt_hit_theta.draw(hsvt_det_theta, "same");
		
		JFrame frame_svt_hit_theta_vs_phi = new JFrame("SVT Hits #theta vs. #phi");
		frame_svt_hit_theta_vs_phi.setSize(1500, 500);
		EmbeddedCanvas can_svt_hit_theta_vs_phi = new EmbeddedCanvas();
		frame_svt_hit_theta_vs_phi.add(can_svt_hit_theta_vs_phi);
		frame_svt_hit_theta_vs_phi.setLocationRelativeTo(null);
		frame_svt_hit_theta_vs_phi.setVisible(true);
		can_svt_hit_theta_vs_phi.divide(1, 1);
		can_svt_hit_theta_vs_phi.cd(0);
		can_svt_hit_theta_vs_phi.setFont("Arial");
		hsvt_det_theta_vs_phi.setTitle("SVT #theta vs #phi");
		hsvt_det_theta_vs_phi.setTitleX("#phi [#degree]");
		hsvt_det_theta_vs_phi.setTitleY("#theta [#degree]");
		can_svt_hit_theta_vs_phi.getPad(0).setTitleFontSize(32);
		can_svt_hit_theta_vs_phi.getPad(0).setAxisTitleFontSize(32);
		can_svt_hit_theta_vs_phi.getPad(0).setAxisLabelFontSize(24);
		can_svt_hit_theta_vs_phi.getPad(0).setStatBoxFontSize(18);
		can_svt_hit_theta_vs_phi.draw(hsvt_det_theta_vs_phi, "same");
		
		JFrame frame_svt_hit_theta_vs_phi_prj = new JFrame("SVT Hits #phi Projections");
		frame_svt_hit_theta_vs_phi_prj.setSize(1000, 500);
		EmbeddedCanvas can_svt_hit_theta_vs_phi_prj = new EmbeddedCanvas();
		frame_svt_hit_theta_vs_phi_prj.add(can_svt_hit_theta_vs_phi_prj);
		frame_svt_hit_theta_vs_phi_prj.setLocationRelativeTo(null);
		frame_svt_hit_theta_vs_phi_prj.setVisible(true);
		can_svt_hit_theta_vs_phi_prj.divide(2, 1);
		can_svt_hit_theta_vs_phi_prj.cd(0);
		can_svt_hit_theta_vs_phi_prj.setFont("Arial");
		hsvt_det_theta_vs_phi_rot.setTitle("Superpositioned SVT #theta vs. #phi");
		hsvt_det_theta_vs_phi_rot.setTitleX("#phi [#degree]");
		hsvt_det_theta_vs_phi_rot.setTitleY("#theta [#degree]");
		can_svt_hit_theta_vs_phi_prj.getPad(0).setTitleFontSize(32);
		can_svt_hit_theta_vs_phi_prj.getPad(0).setAxisTitleFontSize(32);
		can_svt_hit_theta_vs_phi_prj.getPad(0).setAxisLabelFontSize(24);
		can_svt_hit_theta_vs_phi_prj.getPad(0).setStatBoxFontSize(18);
		can_svt_hit_theta_vs_phi_prj.draw(hsvt_det_theta_vs_phi_rot, "same");
		can_svt_hit_theta_vs_phi_prj.cd(1);
		can_svt_hit_theta_vs_phi_prj.setFont("Arial");
		hsvt_det_phi_rot_prj.setTitle("Superpositioned SVT #phi");
		hsvt_det_phi_rot_prj.setTitleX("#phi [#degree]");
		hsvt_det_phi_rot_prj.setTitleY("Counts");
		hsvt_det_phi_rot_prj.setOptStat(10);
		can_svt_hit_theta_vs_phi_prj.getPad(1).setTitleFontSize(32);
		can_svt_hit_theta_vs_phi_prj.getPad(1).setAxisTitleFontSize(32);
		can_svt_hit_theta_vs_phi_prj.getPad(1).setAxisLabelFontSize(24);
		can_svt_hit_theta_vs_phi_prj.getPad(1).setStatBoxFontSize(18);
		can_svt_hit_theta_vs_phi_prj.draw(hsvt_det_phi_rot_prj, "same");
		
		double[] theta_cbnd  = new double[] {40, 59, 63, 66, 68, 70, 72, 74, 77, 80};
		
		JFrame frame_svt_hit_theta_bin = new JFrame("SVT Hits #theta Bins");
		frame_svt_hit_theta_bin.setSize(1000, 1000);
		EmbeddedCanvas can_svt_hit_theta_bin = new EmbeddedCanvas();
		frame_svt_hit_theta_bin.add(can_svt_hit_theta_bin);
		frame_svt_hit_theta_bin.setLocationRelativeTo(null);
		frame_svt_hit_theta_bin.setVisible(true);
		can_svt_hit_theta_bin.divide(3, 3);	
		for(int thbin = 0; thbin < 9; thbin++)
		{
			can_svt_hit_theta_bin.cd(thbin);
			can_svt_hit_theta_bin.setFont("Arial");
			histGroups_svt_det_theta_bin.getItem(thbin).setTitle("SVT " + theta_cbnd[thbin] + "#degree < #theta <= "
																	+ theta_cbnd[thbin+1] + "#degree");
			histGroups_svt_det_theta_bin.getItem(thbin).setTitleX("#theta [#degree]");
			histGroups_svt_det_theta_bin.getItem(thbin).setTitleY("Counts");
			histGroups_svt_det_theta_bin.getItem(thbin).setOptStat(110);
			can_svt_hit_theta_bin.getPad(thbin).setTitleFontSize(32);
			can_svt_hit_theta_bin.getPad(thbin).setAxisTitleFontSize(32);
			can_svt_hit_theta_bin.getPad(thbin).setAxisLabelFontSize(24);
			can_svt_hit_theta_bin.getPad(thbin).setStatBoxFontSize(18);
			can_svt_hit_theta_bin.draw(histGroups_svt_det_theta_bin.getItem(thbin), "same");
		}
		
		JFrame frame_svt_hit_phi_bin = new JFrame("SVT Hits #phi Bins");
		frame_svt_hit_phi_bin.setSize(1000, 1000);
		EmbeddedCanvas can_svt_hit_phi_bin = new EmbeddedCanvas();
		frame_svt_hit_phi_bin.add(can_svt_hit_phi_bin);
		frame_svt_hit_phi_bin.setLocationRelativeTo(null);
		frame_svt_hit_phi_bin.setVisible(true);
		can_svt_hit_phi_bin.divide(2, 2);
		can_svt_hit_phi_bin.cd(0);
		can_svt_hit_phi_bin.setFont("Arial");
		hsvt_det_phi_rot_prj_bin1.setTitle("SVT #phi Bin 1");
		hsvt_det_phi_rot_prj_bin1.setTitleX("#phi [#degree]");
		hsvt_det_phi_rot_prj_bin1.setTitleY("Counts");
		hsvt_det_phi_rot_prj_bin1.setOptStat(110);
		can_svt_hit_phi_bin.getPad(0).setTitleFontSize(32);
		can_svt_hit_phi_bin.getPad(0).setAxisTitleFontSize(32);
		can_svt_hit_phi_bin.getPad(0).setAxisLabelFontSize(24);
		can_svt_hit_phi_bin.getPad(0).setStatBoxFontSize(18);
		can_svt_hit_phi_bin.draw(hsvt_det_phi_rot_prj_bin1, "same");
		can_svt_hit_phi_bin.cd(1);
		can_svt_hit_phi_bin.setFont("Arial");
		hsvt_det_phi_rot_prj_bin2.setTitle("SVT #phi Bin 2");
		hsvt_det_phi_rot_prj_bin2.setTitleX("#phi [#degree]");
		hsvt_det_phi_rot_prj_bin2.setTitleY("Counts");
		hsvt_det_phi_rot_prj_bin2.setOptStat(110);
		can_svt_hit_phi_bin.getPad(1).setTitleFontSize(32);
		can_svt_hit_phi_bin.getPad(1).setAxisTitleFontSize(32);
		can_svt_hit_phi_bin.getPad(1).setAxisLabelFontSize(24);
		can_svt_hit_phi_bin.getPad(1).setStatBoxFontSize(18);
		can_svt_hit_phi_bin.draw(hsvt_det_phi_rot_prj_bin2, "same");
		can_svt_hit_phi_bin.cd(2);
		can_svt_hit_phi_bin.setFont("Arial");
		hsvt_det_phi_rot_prj_bin3.setTitle("SVT #phi Bin 3");
		hsvt_det_phi_rot_prj_bin3.setTitleX("#phi [#degree]");
		hsvt_det_phi_rot_prj_bin3.setTitleY("Counts");
		hsvt_det_phi_rot_prj_bin3.setOptStat(110);
		can_svt_hit_phi_bin.getPad(2).setTitleFontSize(32);
		can_svt_hit_phi_bin.getPad(2).setAxisTitleFontSize(32);
		can_svt_hit_phi_bin.getPad(2).setAxisLabelFontSize(24);
		can_svt_hit_phi_bin.getPad(2).setStatBoxFontSize(18);
		can_svt_hit_phi_bin.draw(hsvt_det_phi_rot_prj_bin3, "same");
		can_svt_hit_phi_bin.cd(3);
		can_svt_hit_phi_bin.setFont("Arial");
		hsvt_det_phi_rot_prj_bin4.setTitle("SVT #phi Bin 4");
		hsvt_det_phi_rot_prj_bin4.setTitleX("#phi [#degree]");
		hsvt_det_phi_rot_prj_bin4.setTitleY("Counts");
		hsvt_det_phi_rot_prj_bin4.setOptStat(110);
		can_svt_hit_phi_bin.getPad(3).setTitleFontSize(32);
		can_svt_hit_phi_bin.getPad(3).setAxisTitleFontSize(32);
		can_svt_hit_phi_bin.getPad(3).setAxisLabelFontSize(24);
		can_svt_hit_phi_bin.getPad(3).setStatBoxFontSize(18);
		can_svt_hit_phi_bin.draw(hsvt_det_phi_rot_prj_bin4, "same");
	}
}