import collections
import gzip
import sys
import itertools
import argparse


class VCF:

    def __init__(self, vcf:str, pop_samples=None):
        self.vcf = vcf
        self.map_index = self.get_index()
        self.samples = self.get_samples()
        self.pop_sample = pop_samples

    @staticmethod
    def get_line(vcf):
        if vcf.endswith('.gz'):
            f = gzip.open(vcf, 'rt')
        else:
            f = open(vcf, 'r')
        for line in f:
            yield line.strip()
        f.close()

    def get_head(self):
        lines = VCF.get_line(self.vcf)
        for line in lines:
            if line.startswith('#CHROM'):
                return line.strip()

    def get_index(self):
        map_index = dict()
        head = self.get_head().split('\t')
        for number, value in enumerate(head):
            map_index[value] = number
        return map_index

    def get_samples(self):
        index = self.get_index()
        return list(index.keys())[9:]

    def __iter__(self):
        for line in VCF.get_line(self.vcf):
            if not line.startswith('#'):
                tags = line.split()
                alt = tags[4]
                if len(alt.split(',')) == 1:
                    if self.pop_sample:
                        yield Record(tags, self.map_index, self.pop_sample)
                    else:
                        yield Record(tags, self.map_index, self.samples)


class FST(VCF):
    INFO = collections.namedtuple('fst', "combine chr_name pos fst")

    def __init__(self, vcf: str, pop_list):
        #继承vcf方法
        super(FST, self).__init__(vcf)
        #执行VCF初始化
        VCF.__init__(self,vcf)
        self.pop_list = pop_list
        self.pop_items = self.pop_dict()
        self.combinations = self.get_combine()

    def get_combine(self):
        return list(itertools.combinations(self.pop_items.keys(), 2))

    def pop_dict(self):
        pop_dict = collections.defaultdict(list)
        with open(self.pop_list) as f:
            for line in f:
                tags = line.strip().split()
                name = tags[0]
                pop = tags[1]
                pop_dict[pop].append(name)
        return pop_dict

    def __iter__(self):
        # all_combine_fst = collections.defaultdict(list)
        for line in VCF.get_line(self.vcf):
            if not line.startswith('#'):
                tags = line.split()
                alt = tags[4]
                if len(alt.split(',')) == 1:
                    for combine in self.combinations:
                        # print(combine)
                        key = '{0}-vs-{1}'.format(*combine)
                        # print(key)
                        group_a = self.pop_items[combine[0]]
                        group_b = self.pop_items[combine[1]]
                        record_a = Record(tags, self.map_index, group_a)
                        record_b = Record(tags, self.map_index, group_b)
                        chr_name = record_a.chrom()
                        pos = record_a.pos()
                        hs_a = record_a.hs()
                        hs_b = record_b.hs()
                        hs = (hs_a + hs_b) / 2
                        alle_cnt_a = record_a.alle_count()
                        alle_cnt_b = record_b.alle_count()
                        all_cnt = sum(alle_cnt_a) + sum(alle_cnt_b)
                        alle_ref_cnt = alle_cnt_a[0] + alle_cnt_b[0]
                        if all_cnt !=0:
                            all_alt_frenquency = alle_ref_cnt / all_cnt
                            ht = 2 * all_alt_frenquency * (1 - all_alt_frenquency)
                            if ht != 0:
                                fst = (ht - hs) / ht
                                yield FST.INFO._make([key, chr_name, pos, fst])
                            else:
                                pass


class PIC:
    INFO = collections.namedtuple('fst', "combine chr_name pos pic")
    def __init__(self, vcf: str, pop_list):
        #继承vcf方法
        super(FST, self).__init__(vcf)
        #执行VCF初始化
        VCF.__init__(self,vcf)
        self.pop_list = pop_list
        self.pop_items = self.pop_dict()

    def pop_dict(self):
        pop_dict = collections.defaultdict(list)
        with open(self.pop_list) as f:
            for line in f:
                tags = line.strip().split()
                name = tags[0]
                pop = tags[1]
                pop_dict[pop].append(name)
        return pop_dict

    def __iter__(self):
        # all_combine_fst = collections.defaultdict(list)
        for line in VCF.get_line(self.vcf):
            if not line.startswith('#'):
                tags = line.split()
                alt = tags[4]
                if len(alt.split(',')) == 1:
                    for pop in self.pop_items.keys():
                        pop_samples = self.pop_items[pop]
                        record = Record(tags, self.map_index, pop_samples)
                        chr_name = record.chrom()
                        pos = record.pos()
                        maf = record.maf()
                        pic = 1 - (pow(maf, 2) + pow((1-maf), 2))
                        yield FST.INFO._make([pop, chr_name, pos, pic])


class Info:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class Record:

    def __init__(self, tags:list, map_index:dict, samples:list):
        self.tags = tags
        self.samples = samples
        self.map_index = map_index

    def samples(self):
        """
        返回样本名称
        :return: list
        """
        print(self.samples)

    @staticmethod
    def get_record(info, sep=';'):
        dict_info = {}
        all_info = info.split(sep)
        for data in all_info:
            map = data.split('=')
            dict_info[map[0]] = dict_info[map[1]]
        return Info(**dict_info)

    def get_format(self):
        """
        获取第八列信息
        :return: info
        """
        Format = self.tags[7]
        return Record.get_record(Format)

    def info(self, sep = ':'):
        info = self.tags[8]
        return info.split(sep)

    def ref(self):
        """
        获取ref的见此类型
        :return: str
        """
        return self.tags[3]

    def alt(self):
        """
        获取突变的基因型
        :return: str
        """
        return self.tags[4]

    def chrom(self):
        """
        获取染色体名称
        :return: str
        """
        return self.tags[0]

    def pos(self):
        """
        获取snp 位置信息
        :return: int
        """
        return self.tags[1]

    def qual(self):
        """
        获取质量值
        :return: float
        """
        return self.tags[5]

    def sample_info(self, sample, sep=':', mistype='./.'):
        """
        获取每个样本的info 信息
        :param sample: 样本名称
        :param sep: 属性之间的分隔符； 例如" "
        :param mistype:
        :return: dict
        """
        sample_index = self.map_index[sample]
        sample_tags = self.tags[sample_index]
        if sample_tags.startswith(mistype) == -1:
            return None
        else:
            sample_info_list = sample_tags.split(sep)
            return dict(zip(self.info(), sample_info_list))

    def get_sample_gt(self, sample, gt='GT'):
        """
        获取个体的基因型
        :param sample:样本名称
        :param gt: 基因型标签表示方式
        :return: str 例如： './.' : '0/1' :'1/1' : '0/0'
        """
        sample_info_dict = self.sample_info(sample)
        if sample_info_dict:
            return sample_info_dict[gt]
        else:
            return None

    def get_sample_ad(self,sample, ad='AD'):
        """
        获取每个体的每个等位基因的深度
        :param sample: 样本名称
        :param ad: 标签AD的表示方式
        :return: str 例如：'1,2'
        """
        sample_info_dict = self.sample_info(sample)
        if sample_info_dict:
            return sample_info_dict[ad]
        else:
            return None

    def get_sample_dp(self, sample, dp='DP'):
        """
        获得每个个体的深度
        :param sample: str,样本名称
        :param dp: 标签深度的标签名称
        :return: int
        """
        sample_info_dict = self.sample_info(sample)
        if sample_info_dict:
            return int(sample_info_dict[dp])
        else:
            return None

    def get_all_gt(self):
        """
        获取每个个体的基因型
        :return: list
        """
        list_gt = list()
        for sample in self.samples:
            gt = self.get_sample_gt(sample)
            if gt:
                list_gt.append(gt)
            else:
                list_gt.append('./.')
        return list_gt

    def get_all_dp(self):
        """
        获取群体中每个个体的深度
        :return: list
        """
        list_dp = list()
        for sample in self.samples:
            dp = self.get_sample_dp(sample)
            if dp:
                list_dp.append(dp)
            else:
                list_dp.append(0)
        return list_dp

    def miss(self, miss_type='./.'):
        """
        计算缺失率
        :param miss_type:
        :return: float
        """
        all_gt = self.get_all_gt()
        miss_count = all_gt.count(miss_type)
        return miss_count/len(all_gt)

    def get_mean_dp(self):
        """
        计算最大深度
        :return: int
        """
        all_dp = self.get_all_dp()
        return max(all_dp)

    def alle_count(self, misstype='./.'):
        """
        返回群体中每个等位基因的频数(过滤掉缺失型)
        :param misstype: 缺失型表示方式
        :return: tuple
        """
        no_miss = ''
        for gt in self.get_all_gt():
            if gt != misstype:
                no_miss += gt
        count_ref = no_miss.count('0')
        cout_alt = no_miss.count('1')
        return (count_ref, cout_alt)

    def maf(self):
        """
        计算最小等位基因频率
        :return: float
        """
        alle_cnt = sum(self.alle_count())
        maf_count = min(self.alle_count())
        if alle_cnt == 0:
            maf = 0
        else:
            maf = maf_count/alle_cnt
        return maf

    def ho(self, misstype='./.',):
        """
        获取观测杂合度
        :param misstype: 缺失型表示方式
        :return: float
        """
        no_miss = []
        ho_cnt = 0
        for gt in self.get_all_gt():
            if gt != misstype:
                no_miss.append(gt)
        for alle in no_miss:
            if alle == '0/1':
                ho_cnt += 1
            elif alle == '1/0':
                ho_cnt += 1
        return ho_cnt/len(no_miss)

    def hs(self):
        """
        计算群体的杂合度
        :return: 2*p*q
        """
        maf = self.maf()
        #2pq
        return 2*(1-maf)*maf

    def get_genotype(self):
        """
        获取每个个体的基因型
        :return: list
        """
        alle_list = list()
        ref = self.ref()
        alt = self.alt()
        fmt = '{}{}'
        for gt in self.get_all_gt():
            if gt == '0/0':
                alle = fmt.format(ref,ref)
            elif (gt == '0/1') or (gt == '1/0'):
                alle = fmt.format(ref,alt)
            elif gt == './.':
                alle = '--'
            elif gt == '1/1':
                alle = fmt.format(alt, alt)
            alle_list.append(alle)
        return alle_list


def get_args():
    parser = argparse.ArgumentParser(description="pase the vcf and calculate some index values")
    parser.add_argument("-v", "--vcf", help="vcf file")
    parser.add_argument("-p", "--pic", help="polymorphism information content")
    parser.add_argument("-f", "--fst", type=str, help="the value of fst")
    parser.add_argument("-list", "--poplist", type=str, help="the file of poplist")
    return parser.parse_args()


def main():
    args = get_args()
    vcf = args.v
    pic = args.pic
    fst = args.f
    pop_list = args.list
    # pop_list = sys.argv[2]
    if fst:
        info = FST(vcf, pop_list)
        # print(info.samples)
        # print(info.map_index)
        for record in info:
            # print(record.sample_info('M'))
            print(record)


main()


































