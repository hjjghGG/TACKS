#include <tna.p4>

#define LOG_MAX_BIN   8
#define EMA_SHIFT     3
#define TOTAL_CELLS   512
#define TOTAL_BINS    (TOTAL_CELLS * LOG_MAX_BIN)
#define PORT_METADATA_SIZE_BYTES 8

typedef bit<32> flow_id_t;
typedef bit<32> ts_t;

header ethernet_t {
    bit<48> dstAddr; bit<48> srcAddr; bit<16> etherType;
}
header ipv4_t {
    bit<4> version; bit<4> ihl; bit<8> dscp_ecn; bit<16> totalLen;
    bit<16> identification; bit<3> flags; bit<13> fragOffset;
    bit<8> ttl; bit<8> protocol; bit<16> hdrChecksum;
    bit<32> srcAddr; bit<32> dstAddr;
}
header tcp_t {
    bit<16> srcPort; bit<16> dstPort; bit<32> seqNo; bit<32> ackNo;
    bit<4> dataOffset; bit<4> res; bit<8> flags; bit<16> window;
    bit<16> checksum; bit<16> urgentPtr;
}
header udp_t {
    bit<16> srcPort; bit<16> dstPort; bit<16> length; bit<16> checksum;
}
header resub_meta_t { bit<32> bin_idx; bit<32> new_center; }

struct ingress_headers_t {
    resub_meta_t resub; ethernet_t ethernet;
    ipv4_t ipv4; tcp_t tcp; udp_t udp;
}
struct port_metadata_t { bit<8> _pad; }

struct ingress_metadata_t {
    flow_id_t flow_id;
    ts_t      now_ts;
    ts_t      last_ts;
    ts_t      gap;
    bit<3>    log_bin_idx;
    bit<32>   bin_idx;
    bit<32>   rd_center;          // old center from register
    bit<32>   new_center;         // EMA result
    bit<32>   diff;
    bit<32>   thresh;
    bit<32>   sv_g;
    bit<16>   sport;
    bit<16>   dport;
};

parser MyIngressParser(packet_in pkt, out ingress_headers_t hdr, out ingress_metadata_t meta,
    out ingress_intrinsic_metadata_t ig_intr_md)
{
    state start {
        meta.bin_idx=0; meta.flow_id=0;
        pkt.extract(ig_intr_md);
        transition select(ig_intr_md.resubmit_flag){
            1w1: parse_resubmit; default: parse_port_metadata;
        }
    }
    state parse_resubmit {
        hdr.resub=pkt.lookahead<resub_meta_t>(); pkt.advance(PORT_METADATA_SIZE_BYTES);
        meta.bin_idx   =hdr.resub.bin_idx;
        meta.new_center=hdr.resub.new_center;
        transition parse_ethernet;
    }
    state parse_port_metadata { pkt.advance(PORT_METADATA_SIZE_BYTES); transition parse_ethernet; }
    state parse_ethernet {
        pkt.extract(hdr.ethernet);
        transition select(hdr.ethernet.etherType){ 0x0800:parse_ipv4; default:accept; }
    }
    state parse_ipv4 {
        pkt.extract(hdr.ipv4);
        transition select(hdr.ipv4.protocol){ 6:parse_tcp; 17:parse_udp; default:accept; }
    }
    state parse_tcp { pkt.extract(hdr.tcp); transition accept; }
    state parse_udp { pkt.extract(hdr.udp); transition accept; }
}

control MyIngress(
    inout ingress_headers_t hdr, inout ingress_metadata_t meta,
    in ingress_intrinsic_metadata_t ig_intr_md,
    in ingress_intrinsic_metadata_from_parser_t ig_prsr_md,
    inout ingress_intrinsic_metadata_for_deparser_t ig_dprsr_md,
    inout ingress_intrinsic_metadata_for_tm_t ig_tm_md)
{
    Register<bit<32>,bit<32>>(TOTAL_CELLS) cell_cnt_reg;
    Register<bit<32>,bit<32>>(TOTAL_CELLS) cell_ts_reg;
    Register<bit<32>,bit<32>>(TOTAL_BINS)  bin_c;
    Register<bit<32>,bit<32>>(TOTAL_BINS)  bin_t;
    Register<bit<32>,bit<32>>(TOTAL_BINS)  bin_m;
    Register<bit<32>,bit<32>>(TOTAL_BINS)  sv_c;
    Register<bit<32>,bit<32>>(TOTAL_BINS)  sv_g;

    // Cell
    RegisterAction<bit<32>,bit<32>,void>(cell_cnt_reg) cc_ra = {
        void apply(inout bit<32> c) { c=c+32w1; }
    };
    RegisterAction<bit<32>,bit<32>,bit<32>>(cell_ts_reg) ct_ra = {
        void apply(inout bit<32> t, out bit<32> p) { p=t; t=meta.now_ts; }
    };

    // Bin center: read old, return old_center (EMA done in MAU)
    RegisterAction<bit<32>,bit<32>,bit<32>>(bin_c) brd_c = {
        void apply(inout bit<32> c, out bit<32> old) { old = c; }
    };
    RegisterAction<bit<32>,bit<32>,void>(bin_c) bwc_ra = {
        void apply(inout bit<32> v) { v = meta.new_center; }
    };

    // Bin total: +1
    RegisterAction<bit<32>,bit<32>,void>(bin_t) bt_ra = {
        void apply(inout bit<32> v) { v = v + 32w1; }
    };

    // Bin match: simple conditional +1 (stability precomputed in MAU)
    RegisterAction<bit<32>,bit<32>,void>(bin_m) bm_ra = {
        void apply(inout bit<32> m) {
            if (meta.sv_g == 32w1) { m = m + 32w1; }
        }
    };

    // Saved center/gap bridge
    RegisterAction<bit<32>,bit<32>,void>(sv_c) sw_c = {
        void apply(inout bit<32> v) { v = meta.rd_center; }
    };
    RegisterAction<bit<32>,bit<32>,bit<32>>(sv_c) sr_c = {
        void apply(inout bit<32> v, out bit<32> r) { r=v; }
    };
    RegisterAction<bit<32>,bit<32>,void>(sv_g) sw_g = {
        void apply(inout bit<32> v) { v = meta.gap; }
    };
    RegisterAction<bit<32>,bit<32>,bit<32>>(sv_g) sr_g = {
        void apply(inout bit<32> v, out bit<32> r) { r=v; }
    };

    Hash<bit<32>>(HashAlgorithm_t.CRC32) flow_hash;
    action set_log_bin(bit<3> v) { meta.log_bin_idx=v; }

    #define T(N,A) table N { actions={A;} default_action=A; size=1; }
    #define A(N,B) action N##_a() B table N { actions={N##_a;} default_action=N##_a; size=1; }

    A(cell_ts,  { meta.last_ts = ct_ra.execute(meta.flow_id); })
    A(cell_cnt, { cc_ra.execute(meta.flow_id); })
    A(bin_s1,   { meta.bin_idx = meta.flow_id << 3; })
    A(bin_s2,   { meta.bin_idx = meta.bin_idx | (bit<32>)meta.log_bin_idx; })
    A(bc_rd,    { meta.rd_center = brd_c.execute(meta.bin_idx); })
    A(ema_s1,   { meta.new_center = meta.rd_center + meta.diff; })      // center + gap3
    A(ema_s2,   { meta.diff = meta.rd_center >> 3; })                      // center>>3
    A(ema_s3,   { meta.new_center = meta.new_center - meta.diff; })        // new_center - center3
    A(bt_inc,   { bt_ra.execute(meta.bin_idx); })
    A(sv_wr_c,  { sw_c.execute(meta.bin_idx); })
    A(sv_wr_g,  { sw_g.execute(meta.bin_idx); })
    A(sv_rd_c,  { meta.rd_center = sr_c.execute(meta.bin_idx); })
    A(sv_rd_g,  { meta.sv_g = sr_g.execute(meta.bin_idx); })
    A(dif_s1,   { meta.diff   = meta.sv_g - meta.rd_center; })
    A(thr_s1,   { meta.thresh = meta.rd_center >> 3; })
    A(stb_s1,   { meta.sv_g = meta.thresh - meta.diff; })            // a = t-d
    A(stb_s2,   { meta.rd_center = meta.diff + meta.thresh; })       // b = d+t
    A(stb_s3,   { meta.sv_g = meta.sv_g >> 31; })                    // sign of a
    A(stb_s4,   { meta.rd_center = meta.rd_center >> 31; })          // sign of b
    A(stb_s5,   { meta.sv_g = meta.sv_g | meta.rd_center; })         // sa|sb
    A(stb_s6,   { meta.sv_g = meta.sv_g ^ 32w1; })                  // stable = not(sa|sb)
    A(mat_s1,   { bm_ra.execute(meta.bin_idx); })
    A(bwc,      { bwc_ra.execute(meta.bin_idx); })

    #undef A
    #undef T

    table gap_bin_table {
        key={meta.gap:ternary;} actions={set_log_bin;}
        const entries={
            (32w0x00001000&&&32w0x00001000):set_log_bin(3w7);
            (32w0x00000400&&&32w0x00000400):set_log_bin(3w6);
            (32w0x00000100&&&32w0x00000100):set_log_bin(3w5);
            (32w0x00000040&&&32w0x00000040):set_log_bin(3w4);
            (32w0x00000010&&&32w0x00000010):set_log_bin(3w3);
            (32w0x00000008&&&32w0x00000008):set_log_bin(3w2);
            (32w0x00000004&&&32w0x00000004):set_log_bin(3w1);
            (32w0x00000000&&&32w0x00000000):set_log_bin(3w0);
        } size=8;
    }

    apply {
        if (ig_intr_md.resubmit_flag == 1w1) {
            sv_rd_c.apply(); sv_rd_g.apply();
            dif_s1.apply();  thr_s1.apply();
            stb_s1.apply(); stb_s2.apply(); stb_s3.apply();
            stb_s4.apply(); stb_s5.apply(); stb_s6.apply();
            mat_s1.apply();
            bwc.apply();
            return;
        }
        meta.now_ts=(bit<32>)ig_prsr_md.global_tstamp[31:0];
        if(hdr.tcp.isValid()){meta.sport=hdr.tcp.srcPort;meta.dport=hdr.tcp.dstPort;}
        else if(hdr.udp.isValid()){meta.sport=hdr.udp.srcPort;meta.dport=hdr.udp.dstPort;}
        else{meta.sport=16w0;meta.dport=16w0;}
        meta.flow_id=flow_hash.get({hdr.ipv4.srcAddr,hdr.ipv4.dstAddr,hdr.ipv4.protocol,meta.sport,meta.dport});
        cell_ts.apply(); cell_cnt.apply();
        meta.gap=meta.now_ts-meta.last_ts; meta.last_ts=(ts_t)meta.last_ts;
        meta.diff = meta.gap >> 3;   // precompute gap3 for EMA
        gap_bin_table.apply();
        bin_s1.apply(); bin_s2.apply();
        bc_rd.apply();
        ema_s1.apply(); ema_s2.apply(); ema_s3.apply();
        bt_inc.apply();           // total+1
        sv_wr_c.apply();          // save old_center
        sv_wr_g.apply();          // save gap
        ig_dprsr_md.resubmit_type=3w1;
    }
}

control MyIngressDeparser(packet_out pkt,inout ingress_headers_t hdr,in ingress_metadata_t meta,
    in ingress_intrinsic_metadata_for_deparser_t ig_dprsr_md)
{
    Resubmit() resubmit;
    apply {
        if(ig_dprsr_md.resubmit_type==3w1){
            resub_meta_t rh; rh.setValid(); rh.bin_idx=meta.bin_idx; rh.new_center=meta.new_center; resubmit.emit(rh);
        }
        pkt.emit(hdr.ethernet);pkt.emit(hdr.ipv4);pkt.emit(hdr.tcp);pkt.emit(hdr.udp);
    }
}

struct egress_headers_t{}; struct egress_metadata_t{};
parser MyEgressParser(packet_in p,out egress_headers_t h,out egress_metadata_t m,
    out egress_intrinsic_metadata_t e){state start{transition accept;}}
control MyEgress(inout egress_headers_t h,inout egress_metadata_t m,
    in egress_intrinsic_metadata_t ei,in egress_intrinsic_metadata_from_parser_t ep,
    inout egress_intrinsic_metadata_for_deparser_t ed,inout egress_intrinsic_metadata_for_output_port_t eo)
{apply{}}
control MyEgressDeparser(packet_out p,inout egress_headers_t h,in egress_metadata_t m,
    in egress_intrinsic_metadata_for_deparser_t ed){apply{}}
Pipeline(MyIngressParser(),MyIngress(),MyIngressDeparser(),
         MyEgressParser(),MyEgress(),MyEgressDeparser()) pipe;
Switch(pipe) main;
