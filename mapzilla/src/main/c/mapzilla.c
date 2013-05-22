#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>
#include <pthread.h>
#include <iostream>
#include <stack>
#include <sparsehash/sparse_hash_map>
#include <stdexcept>
#include <vector>

#include <endian.h>

# if __BYTE_ORDER == __BIG_ENDIAN
/* The host byte order is the same as network byte order,
   so these functions are all just identity.  */
# define ntohll(x)       (x)
# define htonll(x)       (x)
# else
#  if __BYTE_ORDER == __LITTLE_ENDIAN
#   define ntohll(x)     __bswap_64 (x)
#   define htonll(x)     __bswap_64 (x)
#  endif
# endif


using namespace std;
using google::sparse_hash_map;

#define CHUNK_SIZE 25
#define FALSE 0
#define TRUE 1

#define BIG_CONSTANT(x) (x##LLU)

uint64_t MurmurHash64A ( const void * key, int len, uint64_t seed )
{
  const uint64_t m = BIG_CONSTANT(0xc6a4a7935bd1e995);
  const int r = 47;

  uint64_t h = seed ^ (len * m);

  const uint64_t * data = (const uint64_t *)key;
  const uint64_t * end = data + (len/8);

  while(data != end)
  {
    uint64_t k = *data++;

    k *= m;
    k ^= k >> r;
    k *= m;

    h ^= k;
    h *= m;
  }

  const unsigned char * data2 = (const unsigned char*)data;

  switch(len & 7)
  {
  case 7: h ^= uint64_t(data2[6]) << 48;
  case 6: h ^= uint64_t(data2[5]) << 40;
  case 5: h ^= uint64_t(data2[4]) << 32;
  case 4: h ^= uint64_t(data2[3]) << 24;
  case 3: h ^= uint64_t(data2[2]) << 16;
  case 2: h ^= uint64_t(data2[1]) << 8;
  case 1: h ^= uint64_t(data2[0]);
          h *= m;
  };

  h ^= h >> r;
  h *= m;
  h ^= h >> r;

  return h;
}

/*
struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return (s1 == s2) || (s1 && s2 && strncmp(s1, s2, CHUNK_SIZE) == 0);
  }
};

struct my_hash
{
	unsigned long operator()(const char* kmer) const
	{
		unsigned long hash = 0;
		int c;

		for (int i=0; i<CHUNK_SIZE; i++)
		{
			c = kmer[i];
			// hash = hash * 33 ^ c
			hash = ((hash << 5) + hash) ^ c;
		}

		return hash;
	}
};
*/


struct eqstr
{
  bool operator()(long long val1, long long val2) const
  {
	  return val1 == val2;
  }
};

/*
struct my_hash
{
	long long operator()(long long chunk) const
	{
		return chunk;
	}
};
*/

struct my_hash
{
	uint64_t operator()(long long chunk) const
	{
		return MurmurHash64A(&chunk, 8, 97);
		//return chunk;
	}
};

unsigned long hash_func(const char* kmer) {
	unsigned long hash = 0;
	int c;

	for (int i=0; i<CHUNK_SIZE; i++)
	{
		c = kmer[i];
		/* hash = hash * 33 ^ c */
		hash = ((hash << 5) + hash) ^ c;
	}

	return hash;
}

struct location {
	int64_t pos;
//	char ref;
};

struct locations {
	char count;
	char capacity;
	location locs[4];
};

//#define NUM_LOCATIONS 3163842112
#define NUM_LOCATIONS 3200000000UL
//#define NUM_LOCATIONS 200000000
//#define NUM_LOCATIONS 2000000
//#define NUM_LOCATIONS 1000000

int64_t location_idx = 0;
locations* location_pool;

void init_location_pool() {
	printf("sizeof(size_t): %d\n", sizeof(size_t));

	printf("buf size: %lld\n", sizeof(locations) * NUM_LOCATIONS);

	printf("location size: %lld\n", sizeof(locations));

	printf("idx size: %lld\n", sizeof(location_idx));

	fflush(stdout);

	location_pool = (locations*) malloc(sizeof(locations) * NUM_LOCATIONS);
	memset(location_pool, 0, sizeof(locations) * NUM_LOCATIONS);
}

locations* alloc_locations() {
	locations* locs = &(location_pool[location_idx++]);
	locs->count = 0;
	locs->capacity = 4;

	return locs;
}

void add_location(locations* locs, char ref, int64_t pos) {
	if (locs->count < locs->capacity) {
		location* loc = &(locs->locs[locs->count]);
//		loc->ref = ref;
		loc->pos = pos;
		locs->count++;
	}
}

int is_ambiguous(char* chunk) {
	for (int i=0; i<CHUNK_SIZE; i++) {
		if (chunk[i] == 'N') {
			return TRUE;
		}
	}

	return FALSE;
}

sparse_hash_map<long long, locations*, my_hash, eqstr>* chunk_map;
//sparse_hash_map<long long, locations*>* chunk_map;

void print_chunk(const char* prefix, char* str, int idx) {
	char chunk_copy[CHUNK_SIZE+1];
	memset(chunk_copy, 0, CHUNK_SIZE+1);
	memcpy(chunk_copy, str, CHUNK_SIZE);
	printf("%s%s\n", prefix, chunk_copy);
}

locations* find_in_map(long long chunk) {
	locations* locs = NULL;

	sparse_hash_map<long long, locations*, my_hash, eqstr>::const_iterator iter = chunk_map->find(chunk);

	if (iter != chunk_map->end()) {
		locs = (locations*) iter->second;
	}

	return locs;
}

/*
void lookup(char* read, int length) {
	int idx = 0;

	location* chunk_locs[4];
	int chunk_idx = 0;

	while (idx < (length-CHUNK_SIZE)) {
		char* chunk = &(read[idx]);

		char chunk_copy[CHUNK_SIZE+1];
		memset(chunk_copy, 0, CHUNK_SIZE+1);
		memcpy(chunk_copy, chunk, CHUNK_SIZE);

//		printf("Looking up: %s\n", chunk_copy);

//		locations* locs = (*chunk_map)[chunk];
		locations* locs = find_in_map(chunk);

		if (locs != NULL) {
			printf("Found: %d\n", locs->locs[0].pos);
			chunk_locs[chunk_idx++] = &(locs->locs[0]);
//			printf("location: %d\n", locs->locs[0].pos);
		} else {
			printf("Not found: %s\n", chunk_copy);
			chunk_locs[chunk_idx++] = NULL;
//			printf("chunk not found: %d, %s\n", idx, read);
		}
		idx += CHUNK_SIZE;
	}

	int last_pos = -1;

	for (int i=0; i<4; i++) {
		if (chunk_locs[i] != NULL) {

			int pos = chunk_locs[i]->pos - (i*CHUNK_SIZE);

			if (pos == last_pos) {
				printf("READ:\t%d\t%d\t%s\n", chunk_locs[i]->ref, pos+1, read);
				return;
			} else {
				last_pos = pos;
			}
		}
	}

	printf("READ:\t*\t*\t%s\n", read);
}
*/

/*
void inspect_index() {
	printf("Index size: %d\n", chunk_map->size());

	for (sparse_hash_map<long long, locations*, my_hash, eqstr>::const_iterator it = chunk_map->begin();
			 it != chunk_map->end(); ++it) {
//		char* key = (char*) it->first;
		long long key = it->first;
		locations* loc =  it->second;

//		print_chunk("index chunk: ", key, 0);
		printf("key: %lld\n", key);
		printf("index loc1: %d\n", loc->locs[0].pos);
	}
}
*/

void persist_map() {
	printf("Writing index.mz\n");
	FILE* fp = fopen("index.mz", "w");

	if (fp == NULL) {
		printf("Error opening index file");
		exit(-1);
	}

	chunk_map->write_metadata(fp);

	/*
	for (sparse_hash_map<long long, locations*, my_hash, eqstr>::const_iterator it = chunk_map->begin();
			 it != chunk_map->end(); ++it) {
//		char* key = (char*) it->first;
		long long key = it->first;
		locations* loc =  it->second;

//		print_chunk("index chunk: ", key, 0);
		printf("key: %lld\n", key);
		printf("index loc1: %d\n", loc->locs[0].pos);
	}
*/

	printf("sizeof(locations): %d\n", sizeof(locations));

	for (sparse_hash_map<long long, locations*, my_hash, eqstr>::iterator it = chunk_map->begin();
			it != chunk_map->end(); ++it) {

		long long key = it->first;
		locations* value = it->second;

//		printf("key: %lld\n", key);
//		fflush(stdout);
//		printf("v: %x\n", value);
//		fflush(stdout);
//		printf("locs: %x\n", value->locs);
//		fflush(stdout);
//		printf("loc: %d\n", value->locs[0].pos);
//		fflush(stdout);

		fwrite(&key, sizeof(key), 1, fp);
		fwrite(value, sizeof(locations), 1, fp);
	}

	/*
   for (sparse_hash_map<int*, int>::iterator it = ht.begin();
		it != ht.end(); ++it) {
	   const_cast<int*>(it->first) = new int;
	   fread(const_cast<int*>(it->first), sizeof(int), 1, fp);
	   fread(&it->second, sizeof(int), 1, fp);
   }
   */

//	chunk_map-> write_nopointer_data(fp);

	fclose(fp);
	printf("Done writing index.mz\n");
}

void load_map() {
	printf("Reading index.mz\n");
	FILE* fp = fopen("index.mz", "r");

	if (fp == NULL) {
		printf("Error opening index file");
		exit(-1);
	}

	chunk_map->read_metadata(fp);

	printf("Done reading metadata\n");
	fflush(stdout);

	char is_first = TRUE;

	int64_t count = 0;

	for (sparse_hash_map<long long, locations*, my_hash, eqstr>::iterator it = chunk_map->begin();
			it != chunk_map->end(); ++it) {

		if ((count++ % 1000000) == 0) {
			printf("loaded %lld bases from pre-built index\n", count);
		}

//		printf("key val: %lld\n", it->first);
//		printf("key address: %x\n", &(it->first));
//		printf("sizeof first: %d\n", sizeof(it->first));
//		fflush(stdout);

		int ret = fread((void*) &(it->first), sizeof(long long), 1, fp);

//		printf("After fread\n");
//		fflush(stdout);

		if (ret != 1) {
			printf("Error reading key: %d\n", ret);
		}

		locations* locs = alloc_locations();

		ret = fread(locs, sizeof(locations), 1, fp);
		if (ret != 1) {
			printf("Error reading value: %d\n", ret);
		}

		it->second = locs;

//		if (1) {
//			printf("key: %lld\n", it->first);
//			printf("val: %d\n", it->second->locs[0].pos);
//			fflush(stdout);
//			is_first = FALSE;
//		}
	}

	printf("Done loading keys and values\n");
	fflush(stdout);

//	for (sparse_hash_map<int*, int>::iterator it = chunk_map->begin();
//			it != chunk_map->end(); ++it) {
//
//		fread(const_cast<>(it->first), sizeof(int), 1, fp);
//		fread(&it->second, sizeof(int), 1, fp);
//	}


	fclose(fp);
	printf("Done reading index.mz\n");
}

void idx_init() {
	long start = time(NULL);
	printf("Initializing from index...\n");
	chunk_map = new sparse_hash_map<long long, locations*, my_hash, eqstr>();

	printf("Init location pool\n");
	init_location_pool();

	load_map();

	long stop = time(NULL);
	printf("Cache initialized (%ld)\n", (stop-start));
	printf("hash size: %lld\n", chunk_map->size());
	printf("num buckets: %lld\n", chunk_map->bucket_count());
	fflush(stdout);

//	inspect_index();
}

int to2bit(char c) {

	int val;

	switch (c) {
		case 'A':
			val = 0;
			break;
		case 'C':
			val = 1;
			break;
		case 'T':
			val = 2;
			break;
		case 'G':
			val = 3;
			break;
		default:
			printf("Invalid base: %c\n", c);
			exit (-1);
			break;
	}

	return val;
}

long long hash(const char* chunk) {

	long long hash = 0;

	for (int i=0; i<CHUNK_SIZE; i++) {
		char b = chunk[i];
		hash = hash << 2;
		int val = to2bit(b);
		hash += val;
	}

	return hash;
}

long long rc

int init(const char* reference) {

	long start = time(NULL);
	printf("Initializing cache...\n");

	printf("sizetype: %d\n", sizeof(sparse_hash_map<long long, locations*, my_hash, eqstr>::size_type));

	chunk_map = new sparse_hash_map<long long, locations*, my_hash, eqstr>(4000000000UL);
//	chunk_map = new sparse_hash_map<long long, locations*, my_hash, eqstr>(300000000);

	printf("max size: %lld\n", chunk_map->max_size());

//	chunk_map->set_empty_key(NULL);
	printf("Resizing...\n");
	fflush(stdout);
//	chunk_map->resize(155000000);
	printf("Done resizing...\n");
	fflush(stdout);

	FILE *fp = fopen(reference, "r");

	init_location_pool();

	char c = 0;

	char* bases = (char*) malloc(NUM_LOCATIONS);
	memset(bases, 0, NUM_LOCATIONS);
	int64_t idx = 0;

	char is_header = FALSE;

	int num_ambiguous = 0;

	long long count = 0;

	long ssecs = time(NULL);

	while ((c = fgetc(fp)) != EOF) {
//	while (((c = fgetc(fp)) != EOF) && (count++ < 200000)) {

		if ((count++ % 1000000) == 0) {
			long esecs = time(NULL);
			printf("Processed %lld bases.  %ld seconds elapsed.\n", count, (esecs-ssecs));
			printf("hash size: %llu\n", chunk_map->size());
			printf("num buckets: %llu\n", chunk_map->bucket_count());
			fflush(stdout);
			ssecs = time(NULL);
		}

		c = toupper(c);
		if (c == '>') {
			is_header = TRUE;
		} else if (is_header) {
			if (c == '\n') {
				is_header = FALSE;
			}
		} else if (c != '\n') {
			bases[idx++] = c;
			if (idx >= CHUNK_SIZE) {
				char* chunk = &(bases[idx-CHUNK_SIZE]);

				if (!is_ambiguous(chunk)) {
//					hash_func(chunk);

//					locations* locs = (*chunk_map)[chunk];

//					locations* locs = (*chunk_map)[chunk];
					long long chunk_hash = hash(chunk);
					locations* locs = find_in_map(chunk_hash);
					if (locs == NULL) {
						locs = alloc_locations();
						(*chunk_map)[chunk_hash] = locs;
					}

					add_location(locs, '7', idx-CHUNK_SIZE);

//					char chunk_copy[CHUNK_SIZE+1];
//					memset(chunk_copy, 0, CHUNK_SIZE+1);
//					memcpy(chunk_copy, chunk, CHUNK_SIZE);
//					printf("Caching %s\n", chunk_copy);
				} else {
					num_ambiguous++;
				}
			}
		}
	}

	fclose(fp);

	long stop = time(NULL);
	printf("Cache initialized (%ld)\n", (stop-start));
	printf("Num ambiguous positions: %d\n", num_ambiguous);
	printf("hash size: %llu\n", chunk_map->size());
	printf("max hash size: %llu\n", chunk_map->max_size());
	printf("num buckets: %llu\n", chunk_map->bucket_count());
	fflush(stdout);

	persist_map();
//	inspect_index();
}

/*
void proc_stdin() {
	long start = time(NULL);
	printf("Processing start time: %ld\n", start);
	printf("Processing reads\n");
	char read[250];
	memset(read, 0, 250);

	char c = 0;
	int idx = 0;
	while ((c = getchar()) != EOF) {
		if (c != '\n') {
			read[idx++] = c;
		} else {
			read[idx++] = 0;
			lookup(read, idx);
			idx = 0;
		}
	}

	long stop = time(NULL);

	printf("Processing stop(%ld)\n", (stop-start));

}
*/

#define SERVER_PORT 6000

pthread_mutex_t lookup_mutex;

void* do_work(void* t) {

//	long long recvBuf[20000];
	char recvBuff[100000];
	memset(recvBuff, 0, sizeof(recvBuff));

	int* connfd = (int*) t;

	printf("In thread: %d\n", *connfd);
	fflush(stdout);

	char is_done = FALSE;
	while (!is_done) {
		int num_records = 0;
		int n = read(*connfd, &num_records, sizeof(num_records));
		if (n != sizeof(num_records)) {
			printf("Error reading num_records on %d (%d)\n", *connfd, n);
			is_done = TRUE;
		}

//			printf("num_records: %d\n", num_records);
		int num_records2 = ntohl(num_records);
//			printf("num_records2: %d\n", num_records2);
//			fflush(stdout);

		if (num_records2 == -1) {
			is_done = TRUE;
		} else {
			memset(recvBuff, 0, sizeof(recvBuff));
			int size = num_records2 * sizeof(long long);
			int bytes_read = 0;
			char* buf_ptr = recvBuff;
			while (bytes_read < size) {
				n = read(*connfd, buf_ptr, size-bytes_read);
				bytes_read += n;
				buf_ptr = &(recvBuff[bytes_read]);

//					printf("so far: %d\n", bytes_read);
			}

//				printf("finally: %d\n", strlen(recvBuff));

			int64_t positions[100000];
			memset(positions, 0, sizeof(positions));
			//int num_locs = strlen(recvBuff) / CHUNK_SIZE;
			int num_locs = num_records2;

			for (int i=0; i<num_locs; i++) {

				long long* network_chunk_hash = (long long*) &(recvBuff[i * sizeof(long long)]);

				long long chunk_hash = ntohll(*network_chunk_hash);

				locations* locs = find_in_map(chunk_hash);

				if (locs != NULL) {
//					printf("Chunk found: %d\n", locs->locs[0].pos);
					positions[i] = htonll(locs->locs[0].pos);
				} else {
//					printf("Chunk not found %s\n", chunk);
					positions[i] = htonll(-1);
				}
			}

//				int32_t num = htonl(strlen(recvBuff));

//				printf("sending %d, %d\n", num, sizeof(num));

			write(*connfd, positions, sizeof(int64_t) * num_locs);
		}
	}

//		printf("last: %s\n", recvBuff);

	close(*connfd);

	printf("Leaving thread: %d\n", *connfd);
	fflush(stdout);
}

int listen() {
	int listenfd = 0;
	struct sockaddr_in serv_addr;

	int num_threads = 0;
	pthread_t threads[5000];
	int connfd[5000];


	time_t ticks;

	listenfd = socket(AF_INET, SOCK_STREAM, 0);
	memset(&serv_addr, '0', sizeof(serv_addr));

	serv_addr.sin_family = AF_INET;
	serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
	serv_addr.sin_port = htons(SERVER_PORT);

	bind(listenfd, (struct sockaddr*)&serv_addr, sizeof(serv_addr));

	listen(listenfd, 64);

	while(1) {
		printf("Waiting for connection\n");
		fflush(stdout);
		connfd[num_threads] = accept(listenfd, (struct sockaddr*)NULL, NULL);

		printf("Spawning thread: %d\n", connfd[num_threads]);
		fflush(stdout);
		int ret = pthread_create(&(threads[num_threads]), NULL, do_work, &(connfd[num_threads]));
		if (ret != 0) {
			printf("Error creating thread 1: %d\n", ret);
			fflush(stdout);
			exit(-1);
		}

		num_threads++;

		if (num_threads >= 5000) {
			printf("Maximum number of threads reached, exiting\n");
			exit(0);
		}

//		ticks = time(NULL);
//		snprintf(sendBuff, sizeof(sendBuff), "%.24s\r\n", ctime(&ticks));
//		write(connfd, sendBuff, strlen(sendBuff));

//		sleep(1);
	}
}


int main(int argc, char** argv) {

	printf("Starting mapz...\n");

//	map("/datastore/nextgenout2/share/labs/UNCseq/lmose2/mapzilla/chr7.fa");

//	init("/home/lmose/reference/chr7/small7.fa");

//	init("/home/lmose/reference/chr7/small7.fa");
//	proc_stdin();

//	init("/home/lmose/reference/chr7/small7.fa");
//	init("/datastore/nextgenout2/share/labs/UNCseq/lmose2/mapzilla/head7.fa");

	idx_init();

//	init("/datastore/nextgenout2/share/labs/UNCseq/lmose2/mapzilla/chr7.fa");

//	init("/datastore/nextgenout2/share/labs/UNCseq/lmose2/mapzilla/hg19.fa");
	listen();

	printf("mapz done...\n");
}
