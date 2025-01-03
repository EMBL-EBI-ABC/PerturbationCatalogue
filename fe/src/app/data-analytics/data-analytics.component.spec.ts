import { ComponentFixture, TestBed } from '@angular/core/testing';

import { DataAnalyticsComponent } from './data-analytics.component';

describe('DataAnalyticsComponent', () => {
  let component: DataAnalyticsComponent;
  let fixture: ComponentFixture<DataAnalyticsComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [DataAnalyticsComponent]
    })
    .compileComponents();

    fixture = TestBed.createComponent(DataAnalyticsComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
