import React, { useState, useEffect, useRef } from 'react';
import axios from 'axios';
import { 
  Activity, Settings, Terminal, Database, 
  ChevronRight, Play, CheckCircle2, AlertCircle, Loader2,
  Table as TableIcon, Box as CubeIcon, Info, Zap
} from 'lucide-react';
import { clsx } from 'clsx';
import { twMerge } from 'tailwind-merge';

function cn(...inputs) {
  return twMerge(clsx(inputs));
}

// --- Components ---

const SidebarItem = ({ icon: Icon, label, active, onClick }) => (
  <button 
    onClick={onClick}
    className={cn(
      "w-full flex items-center gap-3 px-4 py-3 rounded-lg transition-all duration-200",
      active ? "bg-accent/10 text-accent" : "text-secondary hover:bg-white/5 hover:text-white"
    )}
  >
    <Icon size={20} />
    <span className="font-medium">{label}</span>
  </button>
);

const Card = ({ children, title, icon: Icon, className }) => (
  <div className={cn("glass rounded-xl p-6 flex flex-col gap-4", className)}>
    {title && (
      <div className="flex items-center gap-2 border-b border-white/5 pb-3">
        {Icon && <Icon size={18} className="text-secondary" />}
        <h2 className="font-semibold text-white/90 uppercase tracking-wider text-xs">{title}</h2>
      </div>
    )}
    {children}
  </div>
);

// --- Simple Molstar Loader (Iframe version for stability in this environment) ---
const MolstarIframe = ({ pdbId = '7KOX' }) => (
    <div className="w-full h-full relative rounded-xl overflow-hidden glass border border-white/5">
        <iframe 
            src={`https://www.rcsb.org/3d-view/${pdbId}?preset=cartoon&color=chain`}
            className="w-full h-full border-none bg-black"
            title="Molstar Viewer"
        />
        <div className="absolute bottom-4 right-4 z-10">
            <div className="bg-black/60 backdrop-blur-md px-3 py-1.5 rounded-lg border border-white/10 flex items-center gap-2">
                <CubeIcon size={14} className="text-accent" />
                <span className="text-[10px] font-bold uppercase tracking-widest text-white/80">Interactive 3D View</span>
            </div>
        </div>
    </div>
);

export default function App() {
  const [activeTab, setActiveTab] = useState('pipeline');
  const [logs, setLogs] = useState([]);
  const [steps, setSteps] = useState([]);
  const [running, setRunning] = useState(false);
  const [config, setConfig] = useState({ af_threshold: 0.001, gene_symbol: 'CHRNA7', consequence_filter: 'missense_variant', structure_id: '7kox', species: 'homo_sapiens' });
  const [results, setResults] = useState([]);
  const [currentStep, setCurrentStep] = useState(null);
  
  const logEndRef = useRef(null);
  const ws = useRef(null);

  useEffect(() => {
    fetchConfig();
    fetchSteps();
    loadResults();
  }, []);

  useEffect(() => {
    if (logEndRef.current) {
      logEndRef.current.scrollIntoView({ behavior: 'smooth' });
    }
  }, [logs]);

  const fetchConfig = async () => {
    try {
      const res = await axios.get('http://localhost:8000/config');
      setConfig(res.data);
    } catch (e) {
      console.error("Failed to fetch config", e);
    }
  };

  const fetchSteps = async () => {
    try {
      const res = await axios.get('http://localhost:8000/steps');
      setSteps(res.data);
    } catch (e) {
      console.error("Failed to fetch steps", e);
    }
  };

  const loadResults = async () => {
    try {
      const res = await axios.get(`http://localhost:8000/data/${config.gene_symbol.toLowerCase()}_ranked_variants.csv`);
      setResults(res.data);
    } catch (e) {
      console.error("Failed to load results", e);
    }
  };

  const runPipeline = async () => {
    if (running) return;
    setRunning(true);
    setLogs([]);
    setLogs(prev => [...prev, "Updating parameters and initiating pipeline execution..."]);
    
    try {
      await axios.post('http://localhost:8000/run_pipeline', config);
    } catch (e) {
      setLogs(prev => [...prev, `ERROR: Failed to trigger pipeline: ${e.message}`]);
      setRunning(false);
      return;
    }
    
    ws.current = new WebSocket('ws://localhost:8000/ws/run');
    
    ws.current.onopen = () => {
      ws.current.send(JSON.stringify({ action: 'run_all' }));
    };
    
    ws.current.onmessage = (event) => {
      const data = JSON.parse(event.data);
      if (data.type === 'log') {
        setLogs(prev => [...prev, data.message]);
      } else if (data.type === 'step_start') {
        setCurrentStep(data.index);
      } else if (data.type === 'step_end') {
        // Step finished
      } else if (data.type === 'pipeline_complete') {
        setRunning(false);
        setLogs(prev => [...prev, "Pipeline completed successfully!"]);
        setTimeout(loadResults, 1000);
      } else if (data.type === 'error') {
        setLogs(prev => [...prev, `ERROR: ${data.message}`]);
        setRunning(false);
      }
    };

    ws.current.onclose = () => {
      setRunning(false);
    };
  };

  return (
    <div className="flex w-full h-screen bg-background text-white overflow-hidden">
      {/* Sidebar */}
      <aside className="w-64 border-r border-white/5 flex flex-col p-4 gap-8">
        <div className="flex items-center gap-3 px-2 py-4">
          <div className="w-10 h-10 rounded-lg bg-accent flex items-center justify-center shadow-lg shadow-accent/20">
            <Zap size={24} className="text-white" />
          </div>
          <div>
            <h1 className="font-bold text-lg leading-tight uppercase tracking-tighter">RAREMISS</h1>
            <p className="text-[10px] text-secondary font-bold opacity-60">PIPELINE</p>
          </div>
        </div>

        <nav className="flex-1 flex flex-col gap-1">
          <SidebarItem icon={Activity} label="Dashboard" active={activeTab === 'pipeline'} onClick={() => setActiveTab('pipeline')} />
          <SidebarItem icon={Database} label="Variant Data" active={activeTab === 'data'} onClick={() => setActiveTab('data')} />
          <SidebarItem icon={Settings} label="Configuration" active={activeTab === 'config'} onClick={() => setActiveTab('config')} />
        </nav>

        <div className="p-4 glass rounded-xl flex items-center gap-3">
          <div className={cn("w-2 h-2 rounded-full", running ? "bg-accent animate-pulse" : "bg-success")} />
          <span className="text-xs font-semibold text-secondary uppercase tracking-widest text-[9px]">
            {running ? "Processing" : "Ready"}
          </span>
        </div>
      </aside>

      {/* Main Content */}
      <main className="flex-1 overflow-y-auto bg-gradient-to-br from-transparent to-accent/5 p-8">
        <div className="max-w-7xl mx-auto space-y-8 animate-fade-in">
          
          {activeTab === 'pipeline' && (
            <>
              <div className="flex items-center justify-between">
                <div>
                  <h2 className="text-2xl font-bold tracking-tight">Pipeline Control</h2>
                  <p className="text-secondary text-sm">Orchestrate structural-genomic prioritization routines.</p>
                </div>
                <button 
                  onClick={runPipeline}
                  disabled={running}
                  className={cn(
                    "flex items-center gap-2 px-6 py-3 rounded-lg font-bold transition-all shadow-xl",
                    running ? "bg-white/5 text-secondary cursor-not-allowed" : "bg-accent text-white hover:scale-105 active:scale-95 shadow-accent/20"
                  )}
                >
                  {running ? <Loader2 className="animate-spin" size={20} /> : <Play size={20} />}
                  {running ? "Running Modules..." : "Execute Pipeline"}
                </button>
              </div>

              <div className="grid grid-cols-12 gap-6">
                <div className="col-span-4 flex flex-col gap-6">
                  <Card title="Module Status">
                    <div className="flex flex-col gap-2">
                      {steps.map((step, i) => {
                        const isActive = currentStep === i && running;
                        const isComplete = currentStep > i || (currentStep === i && !running && logs.length > 0);
                        return (
                          <div key={i} className={cn(
                            "flex items-center justify-between p-3 rounded-lg border transition-all",
                            isActive ? "border-accent bg-accent/5" : "border-white/5 bg-white/2 hover:bg-white/5"
                          )}>
                            <div className="flex items-center gap-3">
                              <div className={cn(
                                "w-6 h-6 rounded-full flex items-center justify-center text-[10px] font-bold",
                                isActive ? "bg-accent text-white" : isComplete ? "bg-success/20 text-success" : "bg-white/5 text-secondary"
                              )}>
                                {isComplete ? <CheckCircle2 size={14} /> : i + 1}
                              </div>
                              <span className={cn("text-xs font-medium", isActive ? "text-white" : "text-white/60")}>{step.name}</span>
                            </div>
                            {isActive && <Loader2 size={14} className="animate-spin text-accent" />}
                          </div>
                        );
                      })}
                    </div>
                  </Card>
                </div>

                <div className="col-span-8 flex flex-col gap-6">
                  <Card title="Live Console" icon={Terminal} className="h-[400px] flex flex-col p-4 bg-black/40">
                    <div className="flex-1 overflow-y-auto font-mono text-[11px] leading-relaxed p-2 custom-scrollbar">
                      {logs.length === 0 && <p className="text-secondary opacity-30 italic">Console output will appear here...</p>}
                      {logs.map((log, i) => (
                        <div key={i} className="flex gap-4 border-l border-white/5 pl-4 mb-1">
                          <span className="text-secondary/40 shrink-0 select-none">[{i.toString().padStart(3, '0')}]</span>
                          <span className={cn(
                            log.includes('ERROR') ? "text-danger" : log.includes('SUCCESS') ? "text-success" : "text-white/80"
                          )}>{log}</span>
                        </div>
                      ))}
                      <div ref={logEndRef} />
                    </div>
                  </Card>
                </div>
              </div>
            </>
          )}

          {activeTab === 'data' && (
            <div className="flex flex-col gap-8">
               <div className="flex items-center justify-between">
                <div>
                  <h2 className="text-2xl font-bold tracking-tight text-white">Prioritized Variants</h2>
                  <p className="text-secondary text-sm">Ranked pathogenic candidates with structural evidence.</p>
                </div>
                <button onClick={loadResults} className="text-xs bg-accent/10 border border-accent/20 text-accent px-4 py-2 rounded-lg font-bold hover:bg-accent/20 transition-all">
                    Reload Dataset
                </button>
              </div>

              <div className="grid grid-cols-12 gap-6">
                <div className="col-span-12 lg:col-span-7">
                  <Card icon={TableIcon} className="p-0 overflow-hidden">
                    <div className="max-h-[600px] overflow-y-auto custom-scrollbar">
                      <table className="w-full text-left border-collapse">
                        <thead className="sticky top-0 bg-card/90 backdrop-blur-md z-10">
                          <tr className="border-b border-white/5 text-[10px] uppercase tracking-widest text-secondary font-bold">
                            <th className="px-6 py-4">rsID</th>
                            <th className="px-4 py-4">AA Change</th>
                            <th className="px-4 py-4">Domain</th>
                            <th className="px-4 py-4 text-right">Score</th>
                            <th className="px-6 py-4 text-right">Priority</th>
                          </tr>
                        </thead>
                        <tbody className="text-xs">
                          {results.length === 0 && (
                            <tr>
                              <td colSpan={5} className="px-6 py-24 text-center text-secondary italic opacity-50">
                                No variant records found. Run the prioritization pipeline to populate results.
                              </td>
                            </tr>
                          )}
                          {results.map((row, i) => (
                            <tr key={i} className="border-b border-white/2 hover:bg-accent/5 group transition-all cursor-default">
                              <td className="px-6 py-4 font-bold text-white/90 group-hover:text-accent group-hover:pl-8 transition-all duration-300">{row.rsid}</td>
                              <td className="px-4 py-4 text-white/60">{row.amino_acid_change}</td>
                              <td className="px-4 py-4 text-white/60 italic">{row.domain_region || "-"}</td>
                              <td className="px-4 py-4 text-right font-mono font-bold">{row.priority_score?.toFixed(2)}</td>
                              <td className="px-6 py-4 text-right">
                                <span className={cn(
                                  "px-2 py-1 rounded-md text-[9px] font-black uppercase tracking-widest",
                                  row.priority_category === 'High' ? "bg-danger/20 text-danger border border-danger/30" :
                                  row.priority_category === 'Medium' ? "bg-warning/20 text-warning border border-warning/30" :
                                  "bg-secondary/20 text-secondary border border-secondary/30"
                                )}>
                                  {row.priority_category}
                                </span>
                              </td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  </Card>
                </div>

                <div className="col-span-12 lg:col-span-5 flex flex-col gap-6">
                    <div className="h-[400px]">
                        <MolstarIframe pdbId={config.structure_id} />
                    </div>
                    <Card title="Structural Insights" icon={CubeIcon} className="bg-accent/5 border-accent/10">
                        <div className="space-y-4">
                            <div className="flex gap-3">
                                <div className="p-2 bg-success/20 rounded-lg h-fit"><CheckCircle2 className="text-success" size={16} /></div>
                                <div>
                                    <h4 className="text-sm font-bold text-white/90">Pore Regions</h4>
                                    <p className="text-xs text-secondary leading-relaxed">Variants mapped here are analyzed for their effect on ion conduction and channel gating mechanisms.</p>
                                </div>
                            </div>
                            <div className="flex gap-3">
                                <div className="p-2 bg-accent/20 rounded-lg h-fit"><Info className="text-accent" size={16} /></div>
                                <div>
                                    <h4 className="text-sm font-bold text-white/90">Ligand Binding</h4>
                                    <p className="text-xs text-secondary leading-relaxed">Proximity to nicotine/acetylcholine binding sites is automatically computed in Module 5.</p>
                                </div>
                            </div>
                        </div>
                    </Card>
                </div>
              </div>
            </div>
          )}

          {activeTab === 'config' && (
            <div className="max-w-2xl mx-auto space-y-6">
              <div>
                <h2 className="text-2xl font-bold tracking-tight">System Configuration</h2>
                <p className="text-secondary text-sm">Fine-tune variant filtering and prioritization parameters.</p>
              </div>

              <Card>
                <div className="space-y-6">
                  <div className="grid grid-cols-3 gap-6">
                    <div className="space-y-2">
                      <label className="text-[10px] uppercase font-bold text-secondary tracking-widest">Target Gene</label>
                      <input 
                        type="text" 
                        value={config.gene_symbol}
                        onChange={(e) => setConfig({...config, gene_symbol: e.target.value})}
                        className="w-full bg-black/20 border border-white/10 rounded-lg px-4 py-3 text-sm focus:border-accent outline-none transition-all font-bold"
                      />
                    </div>
                    <div className="space-y-2">
                      <label className="text-[10px] uppercase font-bold text-secondary tracking-widest">Target Species</label>
                      <input 
                        type="text" 
                        value={config.species}
                        onChange={(e) => setConfig({...config, species: e.target.value})}
                        className="w-full bg-black/20 border border-white/10 rounded-lg px-4 py-3 text-sm focus:border-accent outline-none transition-all font-bold"
                      />
                    </div>
                    <div className="space-y-2">
                      <label className="text-[10px] uppercase font-bold text-secondary tracking-widest">Structure ID (PDB)</label>
                      <input 
                        type="text" 
                        value={config.structure_id}
                        onChange={(e) => setConfig({...config, structure_id: e.target.value})}
                        className="w-full bg-black/20 border border-white/10 rounded-lg px-4 py-3 text-sm focus:border-accent outline-none transition-all font-bold uppercase"
                      />
                    </div>
                  </div>
                  <div className="grid grid-cols-2 gap-6">
                    <div className="space-y-2">
                        <label className="text-[10px] uppercase font-bold text-secondary tracking-widest">AF Threshold</label>
                        <input 
                        type="number" 
                        step="0.0001"
                        value={config.af_threshold}
                        onChange={(e) => setConfig({...config, af_threshold: parseFloat(e.target.value)})}
                        className="w-full bg-black/20 border border-white/10 rounded-lg px-4 py-3 text-sm focus:border-accent outline-none transition-all font-mono"
                        />
                    </div>
                    <div className="space-y-2">
                        <label className="text-[10px] uppercase font-bold text-secondary tracking-widest">Mutation Type</label>
                        <select 
                        value={config.consequence_filter}
                        onChange={(e) => setConfig({...config, consequence_filter: e.target.value})}
                        className="w-full bg-black/20 border border-white/10 rounded-lg px-4 py-3 text-sm focus:border-accent outline-none transition-all appearance-none cursor-pointer"
                        >
                            <option value="missense_variant">Missense Variant</option>
                            <option value="stop_gained">Stop Gained</option>
                            <option value="all">All Non-Synonymous</option>
                        </select>
                    </div>
                  </div>

                  <div className="pt-6 border-t border-white/5">
                    <button 
                      onClick={async () => {
                        try {
                          await axios.post('http://localhost:8000/config', config);
                          alert("Configuration saved successfully!");
                        } catch (e) {
                          alert("Failed to save config");
                        }
                      }}
                      className="w-full bg-accent text-white font-bold py-4 rounded-xl hover:brightness-110 active:scale-[0.98] transition-all shadow-lg shadow-accent/20"
                    >
                      Apply Configuration
                    </button>
                  </div>
                </div>
              </Card>

              <div className="bg-accent/5 border border-accent/20 rounded-xl p-6 flex gap-4">
                <Info size={24} className="text-accent shrink-0" />
                <div className="space-y-1">
                  <h4 className="text-sm font-bold text-accent">Framework Notes</h4>
                  <p className="text-xs text-secondary leading-relaxed">
                    Changes to the AF threshold will trigger a fresh query to Ensembl REST during the next Module 1 execution. 
                    Ensure the Target Gene matches the active Structure ID (e.g. {config.structure_id} for {config.gene_symbol}).
                  </p>
                </div>
              </div>
            </div>
          )}
        </div>
      </main>
    </div>
  );
}
