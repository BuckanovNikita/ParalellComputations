#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "mycom.c"
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <pthread.h>

#define MSG_ID   7777
#define MSG_PERM 0600
#define LBUF     100

typedef struct tag_msg_t {
    int n;
    double buf[LBUF];
} msg_t;


int nt, mt;

int msgid;
msg_t msg;

double a = 0;
double b = 1;
int ni = 1000000000;
double sum = 0;

int np, mp;
pthread_t* threads;
double sums[LBUF];

double f1(double x) {
    return 4.0 / (1.0 + x * x);
}

void* myjobt(void* m) {
    int n1, mt = (int) (((long long int) m) % 1024);
    double a1, b1, h1;

    n1 = ni / (nt + np);
    h1 = (b - a) / (nt + np);
    a1 = a + h1 * mt;
    if (mt < (nt + np) - 1) b1 = a1 + h1;
    else b1 = b;

    sums[mt-np] = integrate(f1, a1, b1, n1);

    printf("mt=%d a1=%le b1=%le n1=%d s1=%le\n", mt, a1, b1, n1, sums[mt]);

    return 0;
}


double myjobp(int mp) {
    int n1;
    double a1, b1, h1, s1;

    n1 = ni / (np + nt);
    h1 = (b - a) / (np + nt);
    a1 = a + h1 * mp;
    if (mp < (np+nt) - 1) b1 = a1 + h1;
    else b1 = b;

    s1 = integrate(f1, a1, b1, n1);

    printf("mp=%d a1=%le b1=%le n1=%d s1=%le\n", mp, a1, b1, n1, s1);

    return s1;
}


void NetInit(int np, int* mp);
void NetInit(int np, int* mp) {
    int i;
    pid_t spid = 0;

    if (np > 1) {
        i = 1;
        while (i < np) {
            if (spid > 0 || i == 1) {
                *mp = i;
                spid = fork();
            }
            if (spid == 0) return;
            i++;
        }
    }

    *mp = 0;
    return;
}

int main(int argc, char *argv[]) {
    int i;
    double t;

    if (argc < 3) {
        printf("Usage: %s <process number> <thread number>\n", argv[0]);
        return 1;
    }

    sscanf(argv[1], "%d", &np);
    sscanf(argv[2], "%d", &nt);


    mp = 0;
    mt = 0;

    t = mytime(0);

    if (np + nt < 2)
        sum = integrate(f1, a, b, ni);
    else {
        NetInit(np, &mp);

        if (mp == 0) {
            msgid = msgget(MSG_ID, MSG_PERM | IPC_CREAT);

            if (!(threads = (pthread_t*) malloc(nt * sizeof (pthread_t))))
                myerr("server: not enough memory", 1);

            for (i = 0; i < nt; i++)
                if (pthread_create(threads + i, 0, myjobt, (void*) ((long long int) np + i)))
                    myerr("server: cannot create thread", 2);
        }
        else
            while ((msgid = msgget(MSG_ID, MSG_PERM)) < 0);

        sum = myjobp(mp);

        if (mp > 0) {
            msg.n = 1;
            msg.buf[0] = sum;
            msgsnd(msgid, &msg, sizeof (msg_t), 0);
        } else {
            for (i = 1; i < np; i++) {
                msgrcv(msgid, &msg, sizeof (msg_t), 0, 0);
                sum = sum + msg.buf[0];
            }
            msgctl(msgid, IPC_RMID, (struct msqid_ds *) 0);
        }


        if(mp == 0) {
            for (i = 0; i < nt; i++)
                if (pthread_join(threads[i], 0))
                    myerr("server: cannot wait thread", 3);
                else
                    sum = sum + sums[i];
        free(threads);
        }
    }

    t = mytime(1);

    if (mp == 0)
        printf("time=%lf sum=%22.15le\n", t, sum);

    return 0;
}
