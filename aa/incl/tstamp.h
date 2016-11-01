#ifndef __TSTAMP_H__
#define __TSTAMP_H__

/* ==== HEADER tstamp.h ==== */

/* Useful formatted time routines. */

/* ANSI C, IRIX 5.3, 28. June 1995. Andris Aszodi */

/* ---- PROTOTYPES ---- */

#ifdef __cplusplus
extern "C" {
#endif

/* time_stamp(): returns a pointer to an internal static string
 * that contains a time string in the following format:-
 * "Thu 02-Jun-1994 18:24:23". The string is evaluated at the
 * time of the call and is overwritten between calls.
 */
char *time_stamp(void);

/* start_timer(): starts a timer by saving the current time in an
 * internal variable.
 * Return value: the current calendar time as returned by time().
 */
time_t start_timer(void);

/* stop_timer(): gets the current time and calculates the time
 * elapsed since the first call to start_timer(). Prints a warning
 * if start_timer() hasn't been invoked. Otherwise, returns a ptr to a
 * formatted string that contains the time difference broken down
 * as "26 days 1 hour 3 mins 55 secs".
 */
char *stop_timer(void);

#ifdef __cplusplus
}
#endif

/* ==== END OF HEADER tstamp.h ==== */
#endif
